set(Nornir_SRCS
  itkIRCommon.cxx
  itkIRText.cxx
  itkIRUtils.cxx
  IRTerminator.cxx
  IRThreadStorage.cxx
  IRStdThread.cxx
  IRThreadInterface.cxx
  IRMutex.cxx
  IRMutexInterface.cxx
  IRTransaction.cxx
  IRLog.cxx
  IRThreadPool.cxx
  IRMosaicRefinementCommon.cxx
  IRGridCommon.cxx
  IRFFTCommon.cxx
  IRGridTransform.cxx
  itkRegularStepGradientDescentOptimizer2.cxx
  itkRBFTransform.cxx
  IRAABBox.cxx
  IRTerminator.cxx
  IRStdMutex.cxx
  IRFFT.cxx
)

itk_module_add_library(Nornir ${Nornir_SRCS})
