(ns xl-decoy.digest)

(def aamap {
  \A 71.03712,
  \C 103.00919,
  \D 115.02694,
  \E 129.04259,
  \F 147.06841,
  \G 57.021464,
  \H 137.05891,
  \I 113.08406,
  \K 128.09496,
  \L 113.08406,
  \M 131.04048,
  \N 114.04293,
  \P 97.052764,
  \Q 128.05858,
  \R 156.10111,
  \S 87.032029,
  \T 101.04768,
  \V 99.068414,
  \W 186.07931,
  \Y 163.06333,
  \O 255.1589965,
  \U 144.9595897})

(def WATER (+ (* 2 1.0078250321) 15.99491416))

(defn pepmass [pep mass]
  (let [[fa & oa] pep]
    (if (nil? fa)
      (+ mass WATER)
      (recur oa (+ mass (aamap fa))))))
; (pepmass "AC" 0)

(defn partition-between [f coll]
  (loop [hold [(first coll)]
         rst (rest coll)
         ret []]
    (let [fst (last hold)
          sec (first rst)]
      (cond
        (empty? rst) ret
        (f fst sec) (recur [sec] (rest rst) (conj ret hold))
        :else (recur (conj hold sec) (rest rst) ret)))))

(defn cleaves? [a b spec exc]
  (and (re-matches spec (str a)) (not (re-matches exc (str b)))))

(defn digest [pseq spec exc]
  (let [cle (partition-between (fn [a b] (cleaves? a b spec exc)) pseq)]
    (map clojure.string/join cle)))
; (digest "PEPKEWRERKPAS" #"[KR]" #"[P]")

(defn minpeplen [peps len]
  (filter (fn [item] (>= (count item) len)) peps))

(defn maxpeplen [peps len]
  (filter (fn [item] (<= (count item) len)) peps))

(defn n-missed [peps n]
  (if (<= (count peps) n)
    nil
    (cons (clojure.string/join (apply concat (take (inc n) peps))) (n-missed (next peps) n))))
(defn missed-cleavages [peps n]
  (if (zero? n)
    peps
    (concat (n-missed peps n) (missed-cleavages peps (dec n)))))
