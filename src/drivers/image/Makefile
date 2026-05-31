 ##############################################################################
 # @section DESCRIPTION
 #
 # VIC Image Driver Makefile
 ##############################################################################

# -----------------------------------------------------------------------
# SET ENVIRONMENT-SPECIFIC OPTIONS HERE
# -----------------------------------------------------------------------

# VIC DRIVER PATH (include and src are subdirs of this)
DRIVERPATH = .

# VIC RUN PATH
SHAREDPATH = ../shared_all
SHAREDIMAGEPATH = ../shared_image

# VIC RUN PATH
VICPATH = ../../vic_run

# VIC EXTENTIONS PATH
EXTPATH = ../../extensions
ROUT = rout_rvic

include ${EXTPATH}/${ROUT}/rout.mk

ifndef NC_LIBS
NC_LIBS = $(shell nc-config --libs)
endif

ifndef NC_CFLAGS
NC_CFLAGS = $(shell nc-config --cflags)
endif

GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always --tags)
UNAME_S := $(shell uname -s)
HOSTNAME := $(shell uname -n)


# -------------------------
# LAPACK
# -------------------------
ifeq ($(UNAME_S),Linux)
	LAPACK_INCDIR = /usr/include
	LAPACK_LIBDIR = /usr/lib
	LAPACK_CFLAGS = -I$(LAPACK_INCDIR)
	LAPACK_LIBS   = -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas
else ifeq ($(UNAME_S),Darwin)
    # macOS - Homebrew 安装路径
    LAPACK_INCDIR = /opt/homebrew/opt/lapack/include
    LAPACK_LIBDIR = /opt/homebrew/opt/lapack/lib
	LAPACK_CFLAGS = -I$(LAPACK_INCDIR)
	LAPACK_LIBS   = -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas
endif

ifndef MPICC
MPICC = mpicc
endif
$(info $$MPICC is [${MPICC}])

# Set the log level
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
		   -I ${SHAREDPATH}/include \
		   -I ${SHAREDIMAGEPATH}/include \
		   -I ${EXTPATH}/${ROUT}/include

# Uncomment to include debugging information
CFLAGS  =  ${INCLUDES} ${NC_CFLAGS} ${LAPACK_CFLAGS} -ggdb -O0 -Wall -Wextra -std=c99 \
					 -DLOG_LVL=$(LOG_LVL) \
					 -DGIT_VERSION=\"$(GIT_VERSION)\" \
					 -DUSERNAME=\"$(USER)\" \
					 -DHOSTNAME=\"$(HOSTNAME)\"


LIBRARY = -lm ${NC_LIBS} ${LAPACK_LIBS}

COMPEXE = vic_image
EXT = .exe

# -----------------------------------------------------------------------
# MOST USERS DO NOT NEED TO MODIFY BELOW THIS LINE
# -----------------------------------------------------------------------

HDRS = \
	$(wildcard ${VICPATH}/include/*.h) \
	$(wildcard ${DRIVERPATH}/include/*.h) \
	$(wildcard ${SHAREDPATH}/include/*.h) \
	$(wildcard ${SHAREDIMAGEPATH}/include/*.h) \
	$(wildcard ${EXTPATH}/${ROUT}/include/*.h)

SRCS = \
	$(wildcard ${VICPATH}/src/*.c) \
	$(wildcard ${DRIVERPATH}/src/*.c) \
	$(wildcard ${SHAREDPATH}/src/*.c) \
	$(wildcard ${SHAREDIMAGEPATH}/src/*.c) \
	$(wildcard ${EXTPATH}/${ROUT}/src/*.c)

OBJS = $(SRCS:%.o=%.c)

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
	$(MPICC) -o ${COMPEXE}${EXT} $(OBJS) $(CFLAGS) $(LIBRARY)

# -------------------------------------------------------------
# tags
# so we can find our way around
# -------------------------------------------------------------
tags:	TAGS
TAGS:	$(SRCS) $(HDRS)
	etags $(SRCS) $(HDRS)
clean::
	\rm -f TAGS


# -------------------------------------------------------------
# depend
# -------------------------------------------------------------
depend: .depend
.depend: $(SRCS) $(HDRS)
	$(MPICC) $(CFLAGS) -M $(SRCS) > $@

clean::
	\rm -f .depend
