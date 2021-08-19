test_that("Loading of Platypus VCF files works", {

  example_vcf = system.file(
    "extdata",
    "platypus.vcf.bgz",
    package = "THmisc",
    mustWork = TRUE
  )

  expect_equal(THmisc::indentify_vcf_type(example_vcf), "platypus")
  d1 = expect_silent(THmisc::load_vcf_file(example_vcf, verbose=FALSE))
  expect_s4_class(d1, "CollapsedVCF")

  d2 = expect_output(THmisc::load_vcf_file(example_vcf, verbose=TRUE))
  expect_s4_class(d1, "CollapsedVCF")

})

