(ns xl-decoy.decoy
  (:require [xl-decoy.digest :refer [pepmass digest maxpeplen minpeplen missed-cleavages]]))

(defn twopeps [pep]
  (->> (for [p1 pep p2 pep]
         (concat p1 p2))
       (map clojure.string/join)))

(defn twopepmasses [pep]
  (map (fn [px] (pepmass px 0)) (twopeps pep)))
; (twopepmasses (minpeplen pep 4))

(defn rndprot [pseq]
  (->> (clojure.string/split pseq #"") (shuffle) (clojure.string/join)))

(defn tolcmp [a b]
  (if (< (Math/abs (- a b)) 0.01)
    0
    (if (< a b)
      -1
      1)))

(defn java-bsearch  [xs x]
    (java.util.Collections/binarySearch xs x tolcmp))
 
(defn too-similar [tmass dmass]
  (let [maxct (inc (* (count tmass) 0.0001))
        simcount
        (fn [dm tm ct]
          (if (or (empty? dm) (>= ct maxct))
            ct
            (recur (rest dm) tm (+ ct (if (>= (java-bsearch tm (first dm)) 0) 1 0)))))]
    (< maxct (simcount tmass dmass 0))))

(defn grep 
  ([pattern in] (grep pattern in false))
  ([pattern in invert]
   (let [fun (if invert remove filter)]
     (fun (fn [x] (re-matches pattern (str x))) in))))
;(grep #".*test.*" '("matchingteststring" "nonmatching"))

(defn list-files [dir pattern]
  (let [files (->> dir
                  (clojure.java.io/file)
                  (file-seq)
                  (map str))]
    (grep pattern files)))
;(list-files "resources" #".*txt")

(defn readseq [file]
  (let [oneseq (fn [x] (clojure.string/join
                         (grep #"[^ACDEFGHIKLMNPQRSTUVWY]" (slurp x) true)))]
    (if (= (type file) java.lang.String)
      (list (oneseq file))
      (map oneseq file))))
;(readseq "resources/sequence.txt")

(defn makedecoy [file maxiter]
  (let [pseq (readseq file)
        mkpep (fn [sq] (-> sq
                           (digest #"[KR]" #"[P]")
                           (missed-cleavages 1)
                           (minpeplen 5)
                           (maxpeplen 20)))
        tgt (->> (for [pi (map mkpep pseq)] (twopepmasses pi))
                 (flatten)
                 (sort-by +)
                 (into []))]
    (for [pi pseq]
      (loop [iter maxiter]
        (if (zero? iter)
          false
          (let [dseq (rndprot pi)
                dcy (mkpep dseq)]
            (if (too-similar tgt (into [] (twopepmasses dcy)))
              (recur (dec iter))
              dseq)))))))
