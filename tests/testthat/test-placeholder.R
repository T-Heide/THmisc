test_that("Loading of Platypus VCF files works", {

  example_vcf = system.file(
    "extdata",
    "platypus.vcf.bgz",
    package = "THmisc",
    mustWork = TRUE
  )

  example_bad_vcf = system.file(
    "extdata",
    "platypus_broken.vcf.bgz",
    package = "THmisc",
    mustWork = TRUE
  )

  library(VariantAnnotation)
  expect_equal(THmisc::indentify_vcf_type(example_vcf), "platypus")
  THmisc::load_vcf_file(example_vcf, verbose=FALSE)
  d1 = expect_silent(THmisc::load_vcf_file(example_vcf, verbose=FALSE))
  expect_s4_class(d1, "CollapsedVCF")
  expect_equal(NROW(d1), 2)

  d2 = expect_output(THmisc::load_vcf_file(example_vcf))
  expect_s4_class(d2, "CollapsedVCF")

  expect_error(expect_output(THmisc::load_vcf_file(example_bad_vcf)))
  expect_error(expect_silent(THmisc::load_vcf_file(example_bad_vcf, verbose=FALSE)))
  bd = expect_s4_class(expect_output(THmisc::load_vcf_file(example_bad_vcf, annot=FALSE)), "CollapsedVCF")

  expect_equal(THmisc:::guess_reference_genome(bd), NA)
  expect_equal(THmisc:::guess_reference_genome(d1), 'BSgenome.Hsapiens.UCSC.hg38')


  dannot = expect_s4_class(THmisc::load_vcf_file(example_vcf, annot=TRUE), "CollapsedVCF")
  expect_true("ANNOTATION" %in% names(elementMetadata(rowRanges(dannot))))

})

