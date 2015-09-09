test_qpcrImpute <- function(){
    data(oncogene2013)
    checkEquals(class(oncogene2013)[1],"qPCRset")
    tst <- qpcrImpute(oncogene2013, groupVars=c("sampleType","treatment"))
    checkEquals(class(tst)[1],"qPCRset")
    checkEquals(pData(oncogene2013),pData(tst))
    checkEqualsNumeric(nrow(tst),77)
    checkEqualsNumeric(ncol(tst),24)
    checkEqualsNumeric(exprs(tst)[1,1],24.24607,tolerance=1.0e-4)
    checkEqualsNumeric(exprs(tst)[20,20],28.17372,tolerance=1.0e-4)
}
