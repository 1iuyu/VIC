 ##############################################################################
 # @section DESCRIPTION
 #
 # VIC Classic Driver Makefile
 ##############################################################################

# -----------------------------------------------------------------------
# SET ENVIRONMENT-SPECIFIC OPTIONS HERE
# -----------------------------------------------------------------------

# =============================================================
# Paths
# =============================================================
DRIVERPATH = .
SHAREDPATH = ../shared_all
VICPATH    = ../../vic_run

# -------------------------
# LAPACK
# -------------------------
LAPACK_ROOT = /opt/homebrew/opt/lapack
LAPACK_INCDIR = $(LAPACK_ROOT)/include
LAPACK_LIBDIR = $(LAPACK_ROOT)/lib
LAPACK_CFLAGS = -I$(LAPACK_INCDIR)
LAPACK_LIBS   = -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas

# =============================================================
# Build info
# =============================================================
GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always --tags)
UNAME_S := $(shell uname -s)
HOSTNAME := $(shell uname -n)

# =============================================================
# Include paths
# =============================================================
# To turn off warning statements, set LOG_LVL >= 30
# | Level     | Numeric value    |
# |---------  |---------------   |
# | ERROR     | Always Active    |
# | WARNING   | < 30             |
# | INFO      | < 20             |
# | DEBUG     | < 10             |
LOG_LVL = 30

# set includes
INCLUDES = -I ${DRIVERPATH}/include \
		   -I ${VICPATH}/include \
		   -I ${SHAREDPATH}/include
# Uncomment for normal optimized code flags (fastest run option)
#CFLAGS  = -O3 -Wall -Wno-unused
# LIBRARY = -lm

# =============================================================
# Compile flags
# =============================================================
CFLAGS = \
		${INCLUDES} \
		${NC_CFLAGS} \
		${LAPACK_CFLAGS} \
		-ggdb -O0 \
		-Wall -Wextra \
		-std=c99 \
		-DLOG_LVL=$(LOG_LVL) \
		-DGIT_VERSION=\"$(GIT_VERSION)\" \
		-DUSERNAME=\"$(USER)\" \
		-DHOSTNAME=\"$(HOSTNAME)\"
ifeq (true, ${TRAVIS})
CFLAGS += -rdynamic -Wl,-export-dynamic
endif
# =============================================================
# Link flags
# =============================================================
LDFLAGS = \
	${LAPACK_LIBS} \
	${NC_LIBS} \
	-lm

# =============================================================
# Executable
# =============================================================
COMPEXE = vic_classic
EXT     = .exe

# =============================================================
# Sources & headers
# =============================================================
HDRS = \
	$(wildcard ${VICPATH}/include/*.h) \
	$(wildcard ${DRIVERPATH}/include/*.h) \
	$(wildcard ${SHAREDPATH}/include/*.h) \

SRCS = \
	$(wildcard ${VICPATH}/src/*.c) \
	$(wildcard ${DRIVERPATH}/src/*.c) \
	$(wildcard ${SHAREDPATH}/src/*.c) \
	$(wildcard ${EXTPATH}/${ROUT}/src/*.c)

OBJS = $(SRCS:%.o=%.c)

# =============================================================
# Targets
# =============================================================
all:
	make depend
	make model

default:
	make depend
	make model

full:
	make clean
	make depend
	make tags
	make model

clean::
	\rm -f core log
	\rm -rf ${COMPEXE}${EXT} ${COMPEXE}${EXT}.dSYM

model: $(OBJS)
	$(MPICC) -o ${COMPEXE}${EXT} $(OBJS) $(CFLAGS) $(LDFLAGS)

# =============================================================
# tags
# =============================================================
tags:	TAGS
TAGS:	$(SRCS) $(HDRS)
	etags $(SRCS) $(HDRS)
clean::
	\rm -f TAGS


# =============================================================
# dependencies
# =============================================================
depend: .depend
.depend: $(SRCS) $(HDRS)
	$(MPICC) $(CFLAGS) -M $(SRCS) > $@

clean::
	\rm -f .depend
