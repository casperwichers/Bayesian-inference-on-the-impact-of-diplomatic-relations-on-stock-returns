library(Matrix)

UK = m_plus_list[[1]][1,]
US = m_plus_list[[1]][2,]
China = m_plus_list[[1]][3,]
India = m_plus_list[[1]][4,]
Iran = m_plus_list[[1]][5,]
Australia = m_plus_list[[1]][6,]
Israel = m_plus_list[[1]][7,]
France = m_plus_list[[1]][8,]
Germany = m_plus_list[[1]][9,]
Russia = m_plus_list[[1]][10,]
Brazil = m_plus_list[[1]][11,]
Canada = m_plus_list[[1]][12,]
Turkey = m_plus_list[[1]][13,]
Korea = m_plus_list[[1]][14,]



for(i in 2:length(m_plus_list)){
  UK <- rbind(UK, m_plus_list[[i]][1,])
  US <- rbind(US, m_plus_list[[i]][2,])
  China <- rbind(China, m_plus_list[[i]][3,])
  India <- rbind(India, m_plus_list[[i]][4,])
  Iran <- rbind(Iran, m_plus_list[[i]][5,])
  Australia <- rbind(Australia, m_plus_list[[i]][6,])
  Israel <- rbind(Israel, m_plus_list[[i]][7,])
  France <- rbind(France, m_plus_list[[i]][8,])
  Germany <- rbind(Germany, m_plus_list[[i]][9,])
  Russia <- rbind(Russia, m_plus_list[[i]][10,])
  Brazil <- rbind(Brazil, m_plus_list[[i]][11,])
  Canada <- rbind(Canada, m_plus_list[[i]][12,])
  Turkey <- rbind(Brazil, m_plus_list[[i]][13,])
  Korea <- rbind(Brazil, m_plus_list[[i]][14,])
}

UK = t(UK)
US = t(US)
China = t(China)
India = t(India)
Iran = t(Iran)
Australia = t(Australia)
Israel = t(Israel)
France = t(France)
Germany = t(Germany)
Russia = t(Russia)
Brazil = t(Brazil)
Canada = t(Canada)
Turkey = t(Turkey)
Korea = t(Korea)

rankMatrix(UK)
rankMatrix(US)
rankMatrix(China)
rankMatrix(India)
rankMatrix(Iran)
rankMatrix(Australia)
rankMatrix(Israel)
rankMatrix()

names = c("UK","US")

for(i in names){
  print(rankMatrix(i))
}
