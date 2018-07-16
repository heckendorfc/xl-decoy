(ns xl-decoy.decoy
  (:require [xl-decoy.digest :refer [pepmass digest maxpeplen minpeplen missed-cleavages]]))

(defn twopeps [pep]
  (->> (for [p1 pep p2 pep]
         (concat p1 p2))
       (map clojure.string/join)))

(defn twopepmasses [pep]
  (map #(pepmass % 0) (twopeps pep)))
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
     (fun #(re-matches pattern (str %)) in))))
;(grep #".*test.*" '("matchingteststring" "nonmatching"))

(defn list-files [dir pattern]
  (let [files (->> dir
                  (clojure.java.io/file)
                  (file-seq)
                  (map str))]
    (grep pattern files)))
;(list-files "resources" #".*txt")

(defn readseq [file]
  (let [f-exist #(.exists (clojure.java.io/as-file %))
        safe-slurp #(if (f-exist %) (slurp %) "")
        oneseq #(clojure.string/join
                   (grep #"[^ACDEFGHIKLMNPQRSTUVWY]" (safe-slurp %) true))]
    (if (= (type file) java.lang.String)
      (list (oneseq file))
      (map oneseq file))))
;(readseq "resources/sequence.txt")

(defn decoy-filename [file]
  (map clojure.string/join (map #(concat % ".decoy") file)))

(defn makedecoy 
  ([file maxiter] (makedecoy file maxiter (decoy-filename file)))
  ([file maxiter dfile]
   (let [pseq (readseq file)
         mkpep (fn [sq] (-> sq
                            (digest #"[KR]" #"[P]")
                            (missed-cleavages 1)
                            (minpeplen 5)
                            (maxpeplen 20)))
         tgt (->> (for [pi (map mkpep pseq)] (twopepmasses pi))
                  (flatten)
                  (sort-by +)
                  (into []))
         single-decoy (fn [p-ind]
                        (let [pi (nth pseq p-ind)
                              old-dseq (->> p-ind 
                                            (nth dfile)
                                            (readseq)
                                            (first))]
                          (if (> (count old-dseq) 0) ; decoy file contains a seq already -- recycle
                            old-dseq
                            (loop [iter maxiter]
                              (if (zero? iter)
                                false
                                (let [dseq (rndprot pi)
                                      dcy (mkpep dseq)]
                                  (if (too-similar tgt (into [] (twopepmasses dcy)))
                                    (recur (dec iter))
                                    dseq)))))))]
     (pmap single-decoy (range 0 (count file))))))

(defn savedecoy [file db]
  (if  (< 0  (count file))
    (let [first-file (first file)
          first-file-data (if (false? first-file) "" first-file)]
      (spit  (first file)  (first db))
      (recur  (rest file)  (rest db)))))
