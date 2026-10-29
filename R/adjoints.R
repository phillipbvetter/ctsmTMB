# log determinant ----------------------------------------
# f(X(p)) = log det (X(p))
# df = Tr[X^{-1} dX]
logdet <- RTMB::ADjoint(
  function(x) {
    dim(x) <- rep(sqrt(length(x)), 2)
    log(det(x))
  },
  function(x, y, dy) {
    dim(x) <- rep(sqrt(length(x)), 2)
    t(RTMB::solve(x)) * dy
  },
  name = "logdet")


# Stationary covariance: solve  A P + P A' + Q = 0  for P  (Q = G G')
#
# Vectorised:  (I %x% A + A %x% I) vec(P) = -vec(Q)
#
# Adjoint: for a scalar J with dJ/dP = W, solve the transposed system
#   (I %x% A' + A' %x% I) vec(S) = vec(W)    i.e.  A' S + S A = W
# and then
#   dJ/dQ = -S
#   dJ/dA = -(S P' + S' P)
lyap_matrix <- function(A) {
  I <- diag(nrow(A))
  kronecker(I, A) + kronecker(A, I)
}
solve_lyapunov <- RTMB::ADjoint(
  function(x) {
    n <- sqrt(length(x) / 2)
    A <- x[1:n^2];
    dim(A) <- c(n, n)
    Q <- x[n^2 + 1:n^2]
    RTMB::solve(lyap_matrix(A), -Q)
  },
  function(x, y, dy) {
    n <- sqrt(length(x) / 2)
    A <- x[1:n^2]; dim(A) <- c(n, n)
    P <- y; dim(P) <- c(n, n)
    S <- RTMB::solve(t(lyap_matrix(A)), dy); dim(S) <- c(n, n)
    c(-(S %*% t(P) + t(S) %*% P), -S)
  },
  name = "solve_lyapunov"
)
