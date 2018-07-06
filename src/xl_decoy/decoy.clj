(ns xl-decoy.decoy
  (:require [xl-decoy.digest :refer [pepmass digest minpeplen missed-cleavages]]))

(defn twopeps [pep]
  (->> (for [p1 pep p2 pep]
         (concat p1 p2))
       (map clojure.string/join)))

(defn twopepmasses [pep]
  (map (fn [px] (pepmass px 0)) (twopeps pep)))
; (twopepmasses (minpeplen pep 4))

(defn rndprot [pseq]
  (->> (clojure.string/split pseq #"") (shuffle) (clojure.string/join)))

(defn too-similar [tmass dmass]
  (->> (for [ti tmass, di dmass]
         (< (Math/abs (- ti di)) 0.01))
       (filter identity)
       (count)
       (< (* (count tmass) 0.01))))

(defn makedecoy [file maxiter]
  (let [pseq (clojure.string/join
               (remove (fn [x] (re-matches #"[^ACDEFGHIKLMNPQRSTUVWY]" (str x)))
                       (slurp file)))
        mkpep (fn [sq] (-> sq
                           (digest #"[KR]" #"[P]")
                           (missed-cleavages 3)
                           (minpeplen 4)))
        tgt (twopepmasses (mkpep pseq))]
    (loop [iter maxiter]
      (if (zero? iter)
        false
        (let [dseq (rndprot pseq)
              dcy (mkpep dseq)]
          (if (too-similar tgt (twopepmasses dcy))
            (recur (dec iter))
            dseq))))))
