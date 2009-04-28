
import liblog as log

print log.VERBOSITY
log.prNot(log.VERB_DEBUG, "Hi debug")
log.VERBOSITY +=1

print log.VERBOSITY
log.prNot(log.VERB_DEBUG, "Hi debug")
log.VERBOSITY +=1

print log.VERBOSITY
log.prNot(log.ERROR, "Hi err")
log.VERBOSITY +=1
