library(QAL)
set.seed(2026)
d <- simulate_qal()
x <- paste0("X", 1:10)
f25 <- qal(d, x, "A", "Y", tau = .25)
f75 <- qal(d, x, "A", "Y", tau = .75)
print(f25)
print(f75)
f25$balance
f25$tuning
predict(f25)
