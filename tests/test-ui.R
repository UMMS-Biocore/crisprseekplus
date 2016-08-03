require(shiny)
library(crisprseekplus)
library(shinyjs)
library(testthat)

test_that("trueFalseFunction works as expected", {
    expect_false(trueFalseFunc(2))
    expect_true(trueFalseFunc(1))
    expect_silent(trueFalseFunc(3))
  })

exampleFile <- system.file("extdata","gRNA.fa", package = "GUIDEseq")

test_that("tests fileInputFunction responds correctly to
          a given input file/ receiving no input", {
    expect_null(fileInputFunc(NULL, NULL))
    expect_equal(fileInputFunc(NULL, exampleFile),
                fileInputFunc(exampleFile, NULL))})

test_that("Tests getLogo displays logo", {
    expect_silent(goLogo <- getLogo())
    expect_true(exists("goLogo"))
    expect_equal(goLogo[[1]][[1]], "img")})


