# simulate meristems


# single parameter combo example ------------------------------------------

M  <- 0.5              # genic mutation multiplier
S  <- 0.1              # nonsyn selection
D  <- 0.5              # dominance nonsyn
S2 <- 0                # synonymous selection
S3 <- 0                # non-coding genic selection
D2 <- 0                # dominance synonymous
D3 <- 0                # dominance non-coding

out <- bulk_slim(
  D = D,
  N = 100,
  S = S,
  S2 = S2,
  S3 = S3,
  M = M,
  D2 = D2,
  D3 = D3,
  n_reps = 10,
  slim_file  = "code/meristems.slim",
  fixed_file = "data/slim_out/fixed.txt"
)

nrow(out)


# across parameters -------------------------------------------------------

## parameter grids
Ss   <- c(0, 0.001, 0.01, 0.03, 0.05, 0.1, 0.5, 1)   # nonsyn selection
Ds   <- c(0, 0.1, 0.5, 1)                            # dominance (nonsyn)
Ms   <- seq(from = 0.25, to = 1, length.out = 5)     # genic mutation multiplier
S2s  <- c(0)                                    # synonymous selection

## global-ish
N_pop <- 100
n_reps <- 4000

## build full factorial
params <- CJ(D = Ds,
             S = Ss,
             M = Ms,
             S2 = S2s)

## add things that are fixed across runs
params[, `:=`(
  S3 = 0,
  D2 = D,     # D = D2
  D3 = 0,
  N  = N_pop,
  n_reps  = n_reps,
  slim_file  = "code/meristems.slim",
  fixed_file = "data/slim_out/fixed.txt"
)]

## run them
res_list <- lapply(seq_len(nrow(params)), function(idx) {
  p <- params[idx]

  out <- bulk_slim(
    D  = p$D,
    N  = p$N,
    S  = p$S,
    S2 = p$S2,
    S3 = p$S3,
    M  = p$M,
    D2 = p$D2,
    D3 = p$D3,
    n_reps  = p$n_reps,
    slim_file  = p$slim_file,
    fixed_file = p$fixed_file
  )

  #message(nrow(out))
  message(
    sprintf("Run %d/%d: D=%.3f S=%.3f S2=%.3f M=%.3f -> n=%s (mut/gen ~ %.3f)",
            idx, nrow(params),
            p$D, p$S, p$S2, p$M,
            if (!is.null(out)) nrow(out) else 0,
            if (!is.null(out) && nrow(out) > 0) nrow(out)/p$n_reps/3 else 0)
  )

  # tag output with parameters for later analysis
  if (!is.null(out) && nrow(out) > 0) {
    out[, `:=`(
      D  = p$D,
      S  = p$S,
      S2 = p$S2,
      S3 = p$S3,
      M  = p$M,
      D2 = p$D2,
      D3 = p$D3
    )]
  }

  out
})

all <- rbindlist(res_list, fill = TRUE)
fwrite(all, "tables/SLiM_muts_factorial.csv")


results<-fread("~/Dropbox/Research/developmental_selection/")
