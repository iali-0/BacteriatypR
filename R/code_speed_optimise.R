m <- matrix(1:10,nrow = 5)
v <- c(1,2)
t(m)
m/v

t(t(m)/v) # 10.2 sec

m %*% diag(1/v) # forever

sweep(m,2,v,"/") #10.3sec

m/rep(v, each = nrow(m)) #6.2sec

m/t(replace(t(m),TRUE,v)) # 10.8sec

m/v[col(m)] # 8.9sec
