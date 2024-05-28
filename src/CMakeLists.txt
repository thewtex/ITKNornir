set(Nornir_SRCS
  itkIRCommon.cxx
  itkIRText.cxx
  itkIRUtils.cxx
  IRTerminator.cxx
  IRThreadStorage.cxx
  IRThread.cxx
  IRThreadInterface.cxx
  IRMutex.cxx
  IRMutexInterface.cxx
  IRTransaction.cxx
  IRLog.cxx
  IRThreadPool.cxx
  )

itk_module_add_library(Nornir ${Nornir_SRCS})
