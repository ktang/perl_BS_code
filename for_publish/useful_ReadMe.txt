1_retain_isMeth_with_dep_cutoff_in_all_lib_v0.0.pl

only cytosines with a depth at lease N (usually N=4) in all libraries were retained
for further analysis


| head -1 | perl -F"\t" -lane 'print join("\n", @F)' | less

| perl -F"\t" -lane ' $a = "$F[0]:$F[1]-$F[2]"; $l =$F[2] -  $F[1] + 1;  print join("\t", ( $a, $l,  $F[10], $F[14],$F[18], $F[22],$F[26]  ) ) ' | less

| perl -F"\t" -lane ' $a = "$F[0]:$F[1]-$F[2]"; $l =$F[2] -  $F[1] + 1;  print join("\t", ( $a, $l,  $F[10], $F[14],$F[18], $F[22],$F[26]  ) ) if ($l >=100 ) ' |  perl -lane 'print if($F[-4] / ($F[-3]+ 0.00001)  >= 1.5   and $F[-2] -$F[-1] >= 10 )' |  less