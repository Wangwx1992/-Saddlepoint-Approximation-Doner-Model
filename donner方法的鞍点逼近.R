# --- 1. 参数与函数设置 ---
rho <- 0.9
pi1 <- 0.5
diff_true <- 0.2
n1 <- n2 <- 100
pi2 <- pi1 + diff_true
n_sim <- 5000
# 计算概率
calculate_p <- function(pi_i, rho) {
  p0 <- (1-pi_i)*(rho*pi_i-pi_i+1)
  p1 <- 2*pi_i*(1-pi_i)*(1-rho)
  p2 <- rho*pi_i+(1-rho)*pi_i^2
  return(c(p0 = p0, p1 = p1, p2 = p2))
}

# 模拟生成数据
generate_data <- function(n_i, pi_i, rho) {
  p <- calculate_p(pi_i, rho)
  counts <- rmultinom(1, size = n_i, prob = p)
  return(c(n0 = counts[1], n1 = counts[2], n2 = counts[3]))
}

# 生成观测数据
#set.seed(123)
data1 <- generate_data(n1, pi1, rho)
data2 <- generate_data(n2, pi2, rho)
n01 <- data1["n0"]; n11 <- data1["n1"]; n21 <- data1["n2"]
n02 <- data2["n0"]; n12 <- data2["n1"]; n22 <- data2["n2"]
n1=31;n2=44
# pi1_hat 和 d_hat
pi1_hat <- (n11 + 2 * n21) / (2 * n1)
S0 <- (n01 + n02);S1<- (n11+n12);S2<- (n21+n22)
rho <- (4 * S0*S2-S1^2) / ((2 *S0+S1 )*(S1+2*S2))
d <- (n12 + 2 * n22) / (2 * n2) - (n11 + 2 * n21) / (2 * n1)
n11 =m[2,1];n21=m[3,1]
n12 =m[2,2];n22=m[3,2]
# --- 2. 鞍点逼近主函数 ---
saddlepoint_prob <- function(D) {
 
 
  # 估计概率
  p11 <- 2*pi1_hat*(1-pi1_hat)*(1-rho)
  p21 <- rho*pi1_hat+(1-rho)*pi1_hat^2
  
  p12 <- 2*(pi1_hat+D)*(1-(pi1_hat+D))*(1-rho)
  p22 <- rho*(pi1_hat+D)+(1-rho)*(pi1_hat+D)^2
  
  # 求解鞍点 t，使 K'_Z(t) = 0
  obj <- function(t) {
    A <- p12 * exp(n1 * t) + p22 * exp(2 * n1 * t) + (1 - p12 - p22)
    A1 <- n1 * p12 * exp(n1 * t) + 2 * n1 * p22 * exp(2 * n1 * t)
    B <- p11 * exp(- n2 * t) + p21 * exp(-2  * n2 * t) + (1 - p11 - p21)
    B1 <- - n2 * p11 * exp(- n2 * t) - 2 *  n2 * p21 * exp(-2 *  n2 * t)
    n2 * A1 / A + n1 * B1 / B -2*n1*n2*d
  }
  
  t_hat <- tryCatch(
    uniroot(obj, interval = c(-2, 2), extendInt = "yes")$root,
    error = function(e) NA
  )
  
  if (is.na(t_hat)) return(NA)
  
  t <- t_hat
  # 计算鞍点逼近所需导数和K
  A <- p12 * exp(n1 * t) + p22 * exp(2 * n1 * t) + (1 - p12 - p22)
  A1 <- n1 * p12 * exp(n1 * t) + 2 * n1 * p22 * exp(2 * n1 * t)
  A2 <- n1^2 * p12 * exp(n1 * t) + 4 * n1^2 * p22 * exp(2 * n1 * t)
  A3 <- n1^3 * p12 * exp(n1 * t) + 8 * n1^3 * p22 * exp(2 * n1 * t)
  A4 <- n1^4 * p12 * exp(n1 * t) + 16 * n1^4 * p22 * exp(2 * n1 * t)
  
  B <- p11 * exp(- n2 * t) + p21 * exp(-2 * n2 * t) + (1 - p11 - p21)
  B1 <- - n2 * p11 * exp(- n2 * t) - 2  * n2 * p21 * exp(-2  * n2 * t)
  B2 <-  n2^2 * p11 * exp(- n2 * t) + 4 *  n2^2 * p21 * exp(-2  * n2 * t)
  B3 <- -  n2^3 * p11 * exp(- n2 * t) - 8  * n2^3 * p21 * exp(-2  * n2 * t)
  B4 <-  n2^4 * p11 * exp(- n2 * t) + 16 *  n2^4 * p21 * exp(-2  * n2 * t)
  
  K <- n2 * log(A) + n1 * log(B) - 2*n1*n2*d*t
  K2 <- n2 * (A2 / A - (A1 / A)^2) + n1 * (B2 / B - (B1 / B)^2)
  K3 <- n2 * (A3 / A - 3 * A1 * A2 / A^2 + 2 * (A1 / A)^3) +
    n1 * (B3 / B - 3 * B1 * B2 / B^2 + 2 * (B1 / B)^3)
  K4 <- n2 * (A4 / A - 4 * A1 * A3 / A^2 - 3 * (A2 / A)^2 +
                12 * A1^2 * A2 / A^3 - 6 * (A1 / A)^4) +
    n1 * (B4 / B - 4 * B1 * B3 / B^2 - 3 * (B2 / B)^2 +
            12 * B1^2 * B2 / B^3 - 6 * (B1 / B)^4)
  
  sigma <- sqrt(K2)
  zeta3 <- K3 / (K2^(3/2))
  zeta4 <- K4 / (K2^2)
  
  lambda <- abs(t) * sigma
  gamma_t <- exp(-abs(t)) / (1 - exp(-abs(t)))
  B0 <- lambda * exp(lambda^2 / 2) * (1 - pnorm(lambda))
  B1 <- -lambda * (B0 - dnorm(0))
  B2 <- lambda^2 * (B0 - dnorm(0))
  B3 <- -(lambda^3 * B0 - (lambda^3 - lambda) * dnorm(0))
  B4 <- lambda^4 * B0 - (lambda^4 - lambda^2) * dnorm(0)
  B5 <- -(lambda^5 * B0 - (lambda^5 - lambda^3 + 3 * lambda) * dnorm(0))
  B6 <- lambda^6 * B0 - (lambda^6 - lambda^4 + 3 * lambda^2) * dnorm(0)
  
  MZ <- exp(K)
  term1 <- MZ * exp(-t * d) / (sigma * (1 - exp(-abs(t))))
  term2 <- B0 + (1/sigma) * (1/abs(t) - gamma_t) * B1 + zeta3/6 * sign(t) * B3
  term3 <- (gamma_t * (0.5 - 1/abs(t)) + gamma_t^2) * B2 / sigma^2 +
    (1/sigma) * (1/abs(t) - gamma_t) * zeta3/6 * sign(t) * B4 +
    zeta4/24 * B5 + zeta3^2/72 * B6
  P <- term1 * (term2 + term3)
  return(P)
}

# --- 3. 画图查看 d 对 P 的影响 ---
# 使用新的 d 范围
# 1. 计算 d-grid 下的 saddlepoint 概率值
d_grid <- seq(-1, d + 1, length.out = 1000)
prob_vals <- sapply(d_grid, saddlepoint_prob)

# 2. 画图直观展示
plot(d_grid, prob_vals, type = "l", col = "blue", lwd = 2,
     xlab = expression(d), ylab = " Saddlepoint Approximation Probability",
     #main = "Saddlepoint Approximation vs Delta"
)
abline(h = 0.025, col = "red", lty = 2)
grid()

# 3. 找到 crossing 点（值从 <0.05 到 >0.05 或反之）
cross_index <- which(diff(sign(prob_vals - 0.025)) != 0)

# 4. 判断是否有两个交点
if (length(cross_index) < 2) {
  cat("⚠️ 没找到两个根，可能函数没有交叉两次 0.05。\n")
} else {
  # 提取两个根对应区间
  d_L1 <- d_grid[cross_index[1]]
  d_R1 <- d_grid[cross_index[1] + 1]
  
  d_L2 <- d_grid[cross_index[2]]
  d_R2 <- d_grid[cross_index[2] + 1]
  
  # 逐个求解根
  root1 <- uniroot(function(d) saddlepoint_prob(d) - 0.025,
                   interval = c(d_L1, d_R1))$root
  root2 <- uniroot(function(d) saddlepoint_prob(d) - 0.025,
                   interval = c(d_L2, d_R2))$root
  
  # 输出结果
  cat(sprintf("✅ 找到两个根：\n  delta_lower = %.5f\n  delta_upper = %.5f\n", root1, root2))
  
  # 图上标出根位置
  points(root1, 0.025, col = "green", pch = 19)
  points(root2, 0.025, col = "green", pch = 19)
}



#存储结果的矩阵 ---
results <- matrix(NA, nrow = n_sim, ncol = 4)
colnames(results) <- c("lower", "upper", "length", "covered")

# --- 3. 主循环 ---
for (i in 1:n_sim) {
  # Step 1: 生成模拟数据
  data1 <- generate_data(n1, pi1, rho)
  data2 <- generate_data(n2, pi2, rho)
  
  n01 <- data1["n0"]; n11 <- data1["n1"]; n21 <- data1["n2"]
  n02 <- data2["n0"]; n12 <- data2["n1"]; n22 <- data2["n2"]
  
  # Step 2: 估计 pi1_hat 与 d_hat
  pi1_hat <- (n11 + 2 * n21) / (2 * n1)
  d <- (n12 + 2 * n22) / (2 * n2) - pi1_hat
  
  # Step 3: 构造 delta 网格 & 计算 saddlepoint P(delta)
  d_grid <- seq(-d - 1, d + 1, length.out = 100)
  prob_vals <- sapply(d_grid, function(d) saddlepoint_prob(d))
  
  # Step 4: 找到 crossing 点（P(d) = 0.025）
  cross_index <- which(diff(sign(prob_vals - 0.025)) != 0)
  
  if (length(cross_index) < 2) {
    results[i, ] <- c(NA, NA, NA, NA)  # 无法计算 CI
    next
  }
  
  # 提取 crossing 点附近的区间
  d_L1 <- d_grid[cross_index[1]]
  d_R1 <- d_grid[cross_index[1] + 1]
  d_L2 <- d_grid[cross_index[2]]
  d_R2 <- d_grid[cross_index[2] + 1]
  
  # 求根（精确计算 CI 边界）
  root1 <- tryCatch(
    uniroot(function(d) saddlepoint_prob(d) - 0.025,
            interval = c(d_L1, d_R1))$root,
    error = function(e) NA
  )
  
  root2 <- tryCatch(
    uniroot(function(d) saddlepoint_prob(d) - 0.025,
            interval = c(d_L2, d_R2))$root,
    error = function(e) NA
  )
  
  # 存储结果
  if (is.na(root1) || is.na(root2)) {
    results[i, ] <- c(NA, NA, NA, NA)
  } else {
    ci_length <- root2 - root1
    covered <- (diff_true >= root1) && (diff_true <= root2)
    results[i, ] <- c(root1, root2, ci_length, covered)
  }
  
  # 打印进度（可选）
  if (i %% 100 == 0) cat("已完成", i, "次模拟\n")
}

# --- 4. 计算统计量 ---
# 覆盖率（排除 NA 情况）
coverage <- mean(results[, "covered"], na.rm = TRUE)
# 平均 CI 长度（排除 NA 情况）
avg_length <- mean(results[, "length"], na.rm = TRUE)
# 失败次数（NA 情况）
failures <- sum(is.na(results[, "lower"]))

cat("\n--- 最终结果 ---\n")
cat("覆盖率:", round(coverage, 4), "\n")
cat("平均 CI 长度:", round(avg_length, 4), "\n")
cat("失败次数（无法计算 CI 的情况）:", failures, "/", n_sim, "\n")

