(ns gillespie.steady-state)

(defn- create-stats [state species]
  (-> state
      (assoc :started-recording? true)
      (assoc :stats (into {} (for [specie species]
                               [specie {:W 0, :mu 0, :M 0, :last-modified 0, :history []}])))))

(defn- update-stats [state species]
  (let [{:keys [mixture time]} state
        species (distinct (filter (set species) (map #(nth % 2) (-> state :executed-rxn :cmds))))]
    (update-in state [:stats]
               (fn [state]
                 (reduce (fn [stats specie]
                           (let [{Wn-1 :W, mu_n-1 :mu, Mn-1 :M, lm :last-modified} (stats specie)
                                 wn (- time lm)
                                 Wn (+ Wn-1 wn)
                                 subs (- (mixture specie) mu_n-1)
                                 R (/ (* subs wn) Wn)]
                             (-> stats
                                 (update-in [specie :history] conj {:W Wn-1, :mu mu_n-1, :M Mn-1})
                                 (update-in [specie :history] #(if (> (count %) 10) (subvec % 1) %))
                                 (assoc-in [specie :last-modified] time)
                                 (assoc-in [specie :W] Wn)
                                 (assoc-in [specie :mu] (+ mu_n-1 R))
                                 (assoc-in [specie :M] (+ Mn-1 (* subs Wn-1 R))))))
                         state species)))))

(def *tolerance* 1E-3)
(def *neg-tolerance* (* -1 *tolerance*))

(defn- steady-state? [state]
  (every? (fn [[specie stats]]
            (let [means (map :mu (:history stats))]
              (and (= (count means) 10)
                   (every? #(and (<= % *tolerance*) (>= % *neg-tolerance*))
                           (map (partial apply -) (partition 2 1 means))))))
          (:stats state)))

(defn steady-state [species & {:keys [eq-time eq-steps]}]
  ;; eq-* is the amount of time or steps to simulate before recording the state
  (fn [state]
    (let [update (fn [state]
                   (let [updated-state (update-stats (if (:started-recording? state)
                                                       state
                                                       (create-stats state species)) species)]
                     (if (steady-state? updated-state)
                       (into updated-state {:finished-recording? true
                                            :init-time (:time state)
                                            :init-step (:num-steps state)})
                       updated-state)))
          still-recording? (not (:finished-recording? state))
          need-update? (cond
                         eq-time (and (>= (:time state) eq-time) still-recording?)
                         eq-steps (and (>= (:num-steps state) eq-steps) still-recording?))]
      (if need-update? (update state) state))))

(defn after-steady-state [& {:keys [time num-steps]}]
  (cond
    time #(if (:init-time %) (>= (- (:time %) (:init-time %)) time))
    num-steps #(if (:init-step %) (>= (- (:num-steps %) (:init-step %)) num-steps))))

(comment ;; Usage
(g/simulate [...]
            {...}
            :until (after-steady-state :time 10)
            :callbacks
            (steady-state :eq-time 10)))

