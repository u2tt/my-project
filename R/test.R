library(Rmpfr)
library(parallel)

pi_chudnovsky_parallel <- function(digits, n_cores = detectCores()) {
  precBits <- ceiling(digits * log2(10)) + 64
  n_terms  <- digits %/% 14 + 2
  
  if (n_terms < n_cores) n_cores <- max(1, n_terms)
  
  block_size <- ceiling(n_terms / n_cores)
  starts <- seq(0, n_terms - 1, by = block_size)
  ends   <- pmin(starts + block_size - 1, n_terms - 1)
  n_blocks <- length(starts)
  
  # ---- 串行算出每个块的起始状态 ----
  M <- mpfr(1, precBits)
  L <- mpfr(13591409, precBits)
  X <- mpfr(1, precBits)
  K <- mpfr(6, precBits)
  
  start_states <- vector("list", n_blocks)
  current_i <- 0
  for (b in 1:n_blocks) {
    while (current_i < starts[b]) {
      i_next <- current_i + 1
      M <- M * (K^3 - 16 * K) / (mpfr(i_next, precBits)^3)
      L <- L + 545140134
      X <- X * mpfr("-262537412640768000", precBits)
      K <- K + 12
      current_i <- i_next
    }
    start_states[[b]] <- list(M = M, L = L, X = X, K = K, i = current_i)
  }
  
  # ---- 定义单块计算函数 ----
  compute_block <- function(state, end_i, precBits) {
    M <- state$M; L <- state$L; X <- state$X; K <- state$K; i <- state$i
    S <- mpfr(0, precBits)
    while (i <= end_i) {
      S <- S + (M * L) / X
      i_next <- i + 1
      M <- M * (K^3 - 16 * K) / (mpfr(i_next, precBits)^3)
      L <- L + 545140134
      X <- X * mpfr("-262537412640768000", precBits)
      K <- K + 12
      i <- i_next
    }
    S
  }
  
  # ---- 使用 PSOCK 集群 ----
  cat("启动集群，使用", n_cores, "个核心...\n")
  cl <- makeCluster(n_cores, type = "PSOCK")
  on.exit(stopCluster(cl)) # 确保退出时关闭集群
  
  # 必须在每个 worker 上加载 Rmpfr 包
  clusterEvalQ(cl, library(Rmpfr))
  
  # 导出必要的变量和函数到 worker
  clusterExport(cl, c("compute_block", "precBits", "ends", "start_states"), envir = environment())
  
  # 并行计算
  results <- parLapply(cl, 1:n_blocks, function(b) {
    compute_block(start_states[[b]], ends[b], precBits)
  })
  
  # 汇总
  S_total <- Reduce(`+`, results)
  C <- 426880 * sqrt(mpfr("10005", precBits))
  C / S_total
}

# 示例运行
digits  <- 500000
n_cores <- detectCores(logical = FALSE) # 建议先用物理核心数试试，比如 M2 用 4 或 8
cat("使用", n_cores, "个核心，计算", digits, "位\n")

t0 <- Sys.time()
p  <- pi_chudnovsky_parallel(digits, n_cores)
t1 <- Sys.time()

cat("耗时：", round(as.numeric(t1 - t0, units = "secs"), 2), "秒\n")