
use "C:/Users/Combo/Desktop/Canonical Anlysis STATA oraz R/zdowolenie.dta",clear

codebook rodzina przyjaciele zdrowie sukces 

centile rodzina przyjaciele zdrowie sukces, level(50)

sum rodzina  przyjaciele zdrowie sukces 

sum wiek2011 dochod edukacja kontakty aktywnosc

graph pie, over(stanciv)

ktau rodzina przyjaciele zdrowie sukces, stats(taub p)  

spearman wiek2011 edukacja stanciv dochod zaleznosc kontakty aktywnosc plec zaufanie

xi:canon ( rodzina przyjaciele zdrowie sukces)(dep_wyglad dep_zapal dep_zdrowie dep_sen dep_meczenie i.plec wiek2011 kontakty aktywnosc edukacja i.zaufanie i.zaleznosc i.stanciv dochod), test(1 2 3 4) stdcoef

estat loadings
estat correlations

net install canred

canred 1
canred 2
canred 3
canred 4
