! pressao em Pa
P=1.e5

! calculo do rv
! DEW TEMP em Kelvin 
TD=296.
esln=(23.6837*td-4947.2325)/(td-35.86)
es= exp(esln)
rv=0.622*es/(p-es)

!TEMP em Kelvin
T=299.
ES=610.78*EXP(17.269*(T-273.16)/(T-35.86))
RS=.622*ES/(P-ES)


print*,t,rs,rv,rv/rs*100


end
