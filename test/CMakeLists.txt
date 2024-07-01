itk_module_test()

set(NornirTests
  itkIRRefineGridTest.cxx
  )

CreateTestDriver(Nornir "${Nornir-Test_LIBRARIES}" "${NornirTests}")

itk_add_test(NAME itkIRRefineGridTest
  COMMAND NornirTestDriver
    --compare
    DATA{Baseline/itkIRRefineGridTestOutput.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkIRRefineGridTestOutput.mha
  itkIRRefineGridTest
    ${ITK_TEST_OUTPUT_DIR}/itkIRRefineGridTestOutput.mha
  )
