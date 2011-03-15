## Gillespie's algorithm implementation in <a href="http://clojure.org">Clojure</a>

This Clojure library lets you stochastically simulate chemical reaction systems.
It's goal is to make the interface as simple and natural as I could for chemists.
Also, it's a very concise library of about 150 lines of code.

#### License

The code is released under the <a href="http://www.gnu.org/licenses/lgpl.html">GNU LGPL</a>.

#### Basic tutorial

From the Clojure REPL execute

    (require '[gillespie.core :as g])

This makes the library available under the namespace `g`. You can change the alias as you want.
Then, to simulate for example a reversible Michaelian reaction execute

    (g/simulate ["r1" X0 + E1 <-> X0-E1 at 1, 1
                 "r2" X0-E1 <-> X1-E1 at 1, 1
                 "r3" X1-E1 <-> X1 + E1 at 1, 1]
                {X0 1000
                 E1 20}
                :num-steps 1000)

The code is probably self-explanatory. We are simulation a system of three reactions for 1000 steps.
This returns a sequence with every state of the system during the simulation.
If you are just concerned with the population levels for each molecular specie,
then just call the `g/get-populations` function over the result.
On the other hand, if you want to know which reaction occurred at which time,
call the `g/get-executed-rxns` function over the result.
This will return a map from time values to population levels.

Besides, if you want to simulate a system by a certain amount of internal simulation time instead
of a number of simulation steps, then you can do

    (g/simulate ["r1" X0 + E1 <-> X0-E1 at 1, 1
                 "r2" X0-E1 <-> X1-E1 at 1, 1
                 "r3" X1-E1 <-> X1 + E1 at 1, 1]
                {X0 1000
                 E1 20}
                :time 10) ;; internal time units, which are the same units as for the stochastic kinetic constants

You can also terminate the simulation giving your own predicate function if you know how the simulation
state is represented.
Also, you may include perturbations in your simulation as is shown next

    (g/simulate ["r1" X0 + E1 <-> X0-E1 at 1, 1
                 "r2" X0-E1 <-> X1-E1 at 1, 1
                 "r3" X1-E1 <-> X1 + E1 at 1, 1]
                {X0 1000
                 E1 20}
                :until #(or (>= (:time %) 10) (>= (:num-steps %) 1000))
                :perturbations
                (at 2 time-units add 1 X0)
                (every 20 events del 1 X1)
                (when X0 = 990 change-rate "r1" 10)
                (when "r2" executes inc-rate "r3" 1))

Perturbations are composed by two parts, the predicate part (at the beginning)
and the modification part (at the end).
Posible modifications functions are `add`, `del`, `change-rate`, `inc-rate`, `dec-rate` and `update-rate`.
The last one takes a function you specifies instead of an quantity, unlike the other ones.

But all this wouldn't be so useful if reactions could only be specified in the simulate call,
because two different calls to simulate might commonly share many reactions.
To this end, you can use the `g/rxns` macro, which let you compose a bigger system from smaller sets
of reactions.
So for example

    (let [E1 (g/rxns ["r1" X0 + E1 <-> X0-E1 at 1, 1
                      "r2" X0-E1 <-> X1-E1 at 1, 1
                      "r3" X1-E1 <-> X1 + E1 at 1, 1])
          E2 (g/rxns ["r4" X1 + E2 <-> X1-E2 at 1, 1
                      "r5" X1-E2 <-> X2-E2 at 1, 1
                      "r6" X2-E2 <-> X2 + E2 at 1, 1])]
      (g/simulate (g/rxns E1 E2
                          ["r7" X2 -> at 1])
                  {X0 100000
                   E1 100
                   E2 100}
                  :time 10))

This got better.
But althought we can use now `E1` and `E2` to avoid code repetition when calling `simulate` multiple times,
we still have some code repetition in the definition of `E1` and `E2`.
To rectify this, we can create reactions templates.

    (let [Ej (g/rxns-template [i j]
               "r1" Xi + Ej <-> Xi-Ej at 1, 1
               "r2" Xi-Ej <-> Xj-Ej at 1, 1
               "r3" Xj-Ej <-> Xj + Ej at 1, 1)]
      (g/simulate (g/rxns (Ej 0 1) (Ej 1 2))
                  {X0 100000
                   E1 100
                   E2 100}
                  :time 10))

Better.
But, we know `i` and `j` are correlated, because `j = i+1`.
So let's suppose we want to specify each set of reactions only in terms of `i`,
the subindex of the metabolite that's being catalyzed.
How do we do that?
Well, you have Clojure functions to help you there

    (let [Ej (g/rxns-template [i j]
               Xi + Ej <-> Xi-Ej at 1, 1
               Xi-Ej <-> Xj-Ej at 1, 1
               Xj-Ej <-> Xj + Ej at 1, 1)
          Xi (fn [i] (Ej i (inc i)))]
      (g/simulate (g/rxns (Xi 0) (Xi 1))
                  {X0 100000
                   E1 100
                   E2 100}
                  :time 10))

Perfect.

#### Feedback

If you have any comments, critics or suggestions for improvement,
please let me know by creating an <a href="https://github.com/rhz/gillespie/issues/">issue</a>.

