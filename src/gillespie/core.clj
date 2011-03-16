(ns gillespie.core
  (:require [clojure.contrib.generic.math-functions :as math]
            [clojure.string :as s]))

(defn- choice [m]
  (let [r (rand (apply + (vals m)))]
    (loop [[[k w] & kws] (seq m)
           sum 0]
      (if (or (nil? kws) (< r (+ sum w))) k
          (recur kws (+ sum w))))))

(defn rxn-commands [rxn]
  (:cmds rxn))

(defn execute-command [mixture [cmd n specie]]
  (update-in mixture [specie] (fnil #(cmd % n) 0)))

(defn apply-rxn [rxn mixture]
  (reduce execute-command mixture (rxn-commands rxn)))

(defn get-molecule [state mol]
  (get-in state [:mixture mol] 0))

(defn get-activities [state] ;; TODO this fn takes ~70% of the time. optimize it!
  (into {} (for [{:keys [lhs rate] :as r} (:rxns state)]
             [r (apply * rate (for [[mol freq] lhs]
                                (* (get-molecule state mol) freq)))])))

(defn step [state]
  (let [activities (get-activities state)
        rxn (choice activities)
        total-activity (apply + (vals activities))]
    (if (= total-activity 0)
      (assoc state :deadlock? true)
      (let [dt (/ (math/log (/ 1 (rand))) total-activity)]
        (-> state
            (update-in [:time] (partial + dt))
            (update-in [:num-steps] inc)
            (update-in [:mixture] (partial apply-rxn rxn))
            (assoc :executed-rxn rxn))))))

(defn- get-commands [lhs-counts rhs-counts]
  (for [mol (distinct (concat (keys lhs-counts) (keys rhs-counts)))
        :let [diff (- (or (rhs-counts mol) 0)
                      (or (lhs-counts mol) 0))
              op (cond
                   (> diff 0) +
                   (< diff 0) -)]
        :when op]
    [op (math/abs diff) mol]))

(defn- compile-rxn [[name rxn [rate1 rate2]]]
  (let [[lhs [arrow] rhs] (partition-by #{"<->" "->"} rxn)
        [lhs-counts rhs-counts] (map frequencies [lhs rhs])
        cmds (get-commands lhs-counts rhs-counts)]
    (cond
      (= arrow  "->") [{:name name :lhs lhs-counts :cmds cmds :rate rate1}]
      (= arrow "<->") [{:name name :lhs lhs-counts :cmds cmds :rate rate1}
                       {:name (str name "-op") :lhs rhs-counts :rate rate2
                        :cmds (get-commands rhs-counts lhs-counts)}])))

(defn simulate* [rxns init opts]
  (let [{:keys [num-steps time while until perturbations]} opts
        state {:time 0 :num-steps 0 :mixture init :rxns (mapcat compile-rxn rxns)}
        sim (iterate (fn [state]
                       (reduce #(%2 %1) (step state) perturbations)) state)]
    (cond
      num-steps (take (inc num-steps) sim)
      time (take-while #(<= (:time %) time) sim)
      while (take-while while sim)
      until (take-while (complement until) sim)
      :else sim)))

(defn- encode-rxn [rxn]
  (:rxn (reduce (fn [state a]
                  (if (number? a)
                    (assoc state :n a)
                    {:n 1 :rxn (concat (:rxn state) (repeat (:n state) (str a)))}))
                {:n 1 :rxn []} rxn)))

(defn- or-comb [& preds]
  (fn [x] (some #(% x) preds)))

;; TODO make regular expressions for sequences... it'll make recognizing
;;      patterns in macros much more easy
(defn rxns* [rxn-set-spec]
  (let [splitted (take-nth 2 (partition-by #{'at} rxn-set-spec))
        rates (map (fn [[rxn1 rxn2]]
                     (cond
                       (some #{ '->} rxn1) (take 1 rxn2)
                       (some #{'<->} rxn1) (take 2 rxn2)))
                   (partition 2 1 splitted))
        rxns (cons (first splitted)
                   (drop-last (map #(drop (count %1) %2) rates (rest splitted))))
        names (map (fn [[name]] (if (string? name) name "")) rxns)
        rxns-seq (map (comp encode-rxn (partial remove (or-comb string? #{'* '+}))) rxns)]
    (vec (map vector names (map vec rxns-seq) (map vec rates)))))

(defn transform-if [pred mod]
  (fn [x] (if (pred x) (mod x) x)))

(defmacro rxns [& rxn-set-specs]
  `(apply concat ~(vec (map (transform-if vector? rxns*) rxn-set-specs))))

(defmacro rxns-template [args & rxn-set-specs]
  (let [arg-patterns (map (comp re-pattern str) args)
        args+patterns (vec (map vector args arg-patterns))]
    `(fn ~args (rxns* (for [a# '~rxn-set-specs
                            :let [new-a# (reduce (fn [sym# [arg# argp#]]
                                                   (s/replace sym# argp# (str arg#)))
                                                 (str a#) ~args+patterns)]]
                        (cond
                          (= a# 'at) 'at
                          (symbol? a#) (symbol new-a#)
                          (string? a#) new-a#
                          (number? a#) a#))))))

(defn- make-greater-than-zero [x]
  (if (< x 0) 0 x))

(defn update-molecule [state specie f]
  (update-in state [:mixture specie] (comp make-greater-than-zero (fnil f 0))))

(defn update-rate [state rxn-name f]
  (update-in state [:rxns]
             (partial map #(if (= (:name %) rxn-name) (update-in % [:rate] f) %))))

(defn- make-modification [[op qty-or-rxn specie-or-qty-or-fn]]
  (case op
    'add #(update-molecule % (str specie-or-qty-or-fn) (partial + qty-or-rxn))
    'del #(update-molecule % (str specie-or-qty-or-fn) (partial - qty-or-rxn))
    'update-rate #(update-rate % qty-or-rxn specie-or-qty-or-fn)
    'change-rate #(update-rate % qty-or-rxn (partial + specie-or-qty-or-fn))
    'inc-rate #(update-rate % qty-or-rxn (partial + specie-or-qty-or-fn))
    'dec-rate #(update-rate % qty-or-rxn (partial - specie-or-qty-or-fn))))

(defn- get-qty-if-molecule [x state]
  (let [{:keys [mixture]} state]
    (if (symbol? x) (mixture (str x)) x)))

(defn- make-pred [op pred-part]
  (case op
    'at (let [[wait-amount wait-type] pred-part]
          #(= wait-amount ((case wait-type
                             'time-units :time
                             'events :num-steps) %)))
    'every (let [[wait-amount wait-type] pred-part]
             #(zero? (mod ((case wait-type
                             'time-units :time
                             'events :num-steps) %) wait-amount)))
    'when (cond
            (= (count pred-part) 3) (let [[a op b] pred-part]
                                      #((resolve op)
                                        (get-qty-if-molecule a %)
                                        (get-qty-if-molecule b %)))
            (= (last pred-part) 'executes) #(= (-> % :executed-rxn :name)
                                               (second pred-part)))))

(defn parse-perturbations [perturbations]
  (for [p perturbations
        :let [mod-part (take-last 3 p)
              pred-part (drop-last 3 (rest p))]]
    (transform-if (make-pred (first p) pred-part)
                  (make-modification mod-part))))

(defn- parse-opts [opts]
  (let [groups (map (partial apply concat) (partition 2 (partition-by keyword? opts)))
        [[_ & perturbations]] (filter (comp #{:perturbations} first) groups)]
    ;; delay parse-perturbations until runtime, because it returns fn objects
    (into {:perturbations `(parse-perturbations ~perturbations)}
          (map vec (remove (comp #{:perturbations} first) groups)))))

(defn replace-keys [f m]
  (into {} (for [[k v] m] [(f k) v])))

(defmacro simulate [rxn-set-spec init & opts]
  `(simulate* ~((transform-if vector? rxns*) rxn-set-spec)
              ~(replace-keys str init) ~(parse-opts opts)))

(defn get-populations [simulation]
  (zipmap (map :time simulation) (map :mixture simulation)))

(defn get-executed-rxns [simulation]
  (into {} (for [{:keys [time executed-rxn]} (rest simulation)]
             [time (:name executed-rxn)])))

