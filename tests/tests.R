library(testthat);
context("EditR ");
source("global.R");

test_that(desc = "Default trimming on example sequencing file", code = {
  
  filename <- "tests/seqfiles/example.ab1"
  guideseq <- "CACTGGAATGACACACGCCC"
  p.val.cutoff = 0.01
  default.trim = TRUE
  is.reverse = FALSE
  
  results <- CalcEditR(filename, guideseq, p.val.cutoff, default.trim, is.reverse)
  expect_that(object=results$trim$trim5, condition=equals(20))
})

# negative tests

test_that(desc = "Error is thrown when the reads are noisy", code = {
  
  filename <- "tests/seqfiles/noisy_reads.ab1"
  guideseq <- "TGAGGCCCAAAGAAGTTTTC"
  p.val.cutoff = 0.01
  default.trim = FALSE
  is.reverse = TRUE
  
  err <- expect_error(CalcEditR(filename, guideseq, p.val.cutoff, default.trim, is.reverse));
})


test_that(desc = "Set reverse complement when it is not supposed to be", code = {
  
  filename <- "tests/seqfiles/example.ab1"
  guideseq <- "CACTGGAATGACACACGCCC"
  p.val.cutoff = 0.01
  default.trim = TRUE
  is.reverse = TRUE
  
  expect_error(CalcEditR(filename, guideseq, p.val.cutoff, default.trim, is.reverse));
  
})

test_that(desc = "Use wrong gRNA", code = {
  
  filename <- "tests/seqfiles/example.ab1"
  guideseq <- "TGAGGCCCAAAGAAGTTTTC"
  p.val.cutoff = 0.01
  default.trim = TRUE
  is.reverse = FALSE
  
  expect_error(CalcEditR(filename, guideseq, p.val.cutoff, default.trim, is.reverse));
  
})


# Test automatic trimming

test_that(desc = "Trim 3' end", code = {
  
  filename <- "tests/seqfiles/noisy_past_230.ab1"
  guideseq <- "GTACAACTCGGCTCCTGTTG"
  p.val.cutoff = 0.01
  default.trim = FALSE
  is.reverse = TRUE
  
  results <- CalcEditR(filename, guideseq, p.val.cutoff, default.trim, is.reverse)
  expect_lt(results$trim$trim3, 232)
  expect_lt(object=results$trim$trim5, 30)
})

test_that(desc = "Trim 3' end with low noise", code = {
  
  filename <- "tests/seqfiles/low_noise.ab1"
  guideseq <- "ATAAATACTTAGTTATGCTG"
  p.val.cutoff = 0.01
  default.trim = FALSE
  is.reverse = TRUE
  
  results <- CalcEditR(filename, guideseq, p.val.cutoff, default.trim, is.reverse)
  expect_lt(results$trim$trim3, 470)
  expect_lt(object=results$trim$trim5, 30)
})

test_that(desc = "Clean read", code = {
  
  filename <- "tests/seqfiles/clean_read.ab1"
  guideseq <- "TGAAGAGCTTCCATCTGATT"
  p.val.cutoff = 0.01
  default.trim = FALSE
  is.reverse = FALSE
  
  results <- CalcEditR(filename, guideseq, p.val.cutoff, default.trim, is.reverse)
  expect_gt(results$trim$trim3, 400)
  expect_lt(object=results$trim$trim5, 30)
})

test_that(desc = "Short read", code = {
  
  filename <- "tests/seqfiles/short_read.ab1"
  guideseq <- "TACATCCCCCCCTCTTGCAA"
  p.val.cutoff = 0.01
  default.trim = FALSE
  is.reverse = FALSE
  
  results <- CalcEditR(filename, guideseq, p.val.cutoff, default.trim, is.reverse)
  expect_lt(results$trim$trim3, 240)
  expect_gt(results$trim$trim3, 220)
  expect_lt(object=results$trim$trim5, 30)
})
