test_that("EPICC annotations work", {
  exp_annot = THmisc:::epicc_wgs_annots
  expect_identical(annotation_from_barcode(exp_annot$sample_barcode), exp_annot)
})

