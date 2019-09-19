.PHONY: clean static deps .prebuild

CWD:=$(shell pwd)
SRC_DIR:=$(CWD)/src
BIN_DIR:=$(CWD)/bin
OBJ_DIR:=$(CWD)/obj
LIB_DIR:=$(CWD)/lib
DEP_DIR:=$(CWD)/deps
INC_DIR:=$(CWD)/include

EXE:=acgtn
all: $(BIN_DIR)/$(EXE)

CXX = g++
CXX_FLAGS += -std=c++11 -O3
CXX_FLAGS += -DNDEBUG -DBITS_INDEX_32
STATIC_FLAGS=-static -static-libstdc++ -static-libgcc
LD_INCLUDE_FLAGS:=-I$(INC_DIR) -I$(SRC_DIR)
LD_LIB_FLAGS:= -L$(LIB_DIR) -lrocksdb -lpfor -lzstd -lsnappy -llz4 -llzma -lbz2 -lz -lpthread -Wl,-rpath,$(LIB_DIR)

CXX_OBJS = $(filter-out $(OBJ_DIR)/main.o,$(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(wildcard $(SRC_DIR)/*.cpp)))

XZ_DIR:=$(DEP_DIR)/xz-5.2.4
ZLIB_DIR:=$(DEP_DIR)/zlib-1.2.11
BZIP2_DIR:=$(DEP_DIR)/bzip2-1.0.6
LZ4_DIR:=$(DEP_DIR)/lz4-1.8.0
SNAPPY_DIR:=$(DEP_DIR)/snappy-1.1.4
ZSTD_DIR:=$(DEP_DIR)/zstd-1.4.0
ROCKSDB_DIR:=$(DEP_DIR)/rocksdb-6.1.2
PFOR_DIR:=$(DEP_DIR)/TurboPFor-master
ROBINMAP_DIR:=$(DEP_DIR)/robin-map-0.6.1

LIB_DEPS =
LIB_DEPS += $(LIB_DIR)/liblzma.a
LIB_DEPS += $(LIB_DIR)/libz.a
LIB_DEPS += $(LIB_DIR)/libbz2.a
LIB_DEPS += $(LIB_DIR)/liblz4.a
LIB_DEPS += $(LIB_DIR)/libsnappy.a
LIB_DEPS += $(LIB_DIR)/libzstd.a
LIB_DEPS += $(LIB_DIR)/librocksdb.a
LIB_DEPS += $(LIB_DIR)/libpfor.a
LIB_DEPS += copytsl

CXX_ENV = export LD_LIBRARY_PATH=$(LIB_DIR):$(LD_LIBRARY_PATH); export LD_INCLUDE_PATH=$(INC_DIR):$(LD_INCLUDE_PATH);

$(BIN_DIR)/acgtn: $(CXX_OBJS) $(LIB_DEPS)
	$(CXX_ENV)$(CXX) $(CXX_FLAGS) -o $(BIN_DIR)/$(EXE) $(SRC_DIR)/main.cpp $(CXX_OBJS) $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS)

static: $(CXX_OBJS) $(LIB_DEPS)
	$(CXX_ENV)$(CXX) $(CXX_FLAGS) $(STATIC_FLAGS) -o $(BIN_DIR)/$(EXE) $(SRC_DIR)/main.cpp $(CXX_OBJS) $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS)

deps: $(LIB_DEPS)

$(LIB_DIR)/libz.a:
	$(CXX_ENV)cd $(ZLIB_DIR) && CFLAGS='-fPIC ${EXTRA_CFLAGS}' LDFLAGS='${EXTRA_LDFLAGS}' ./configure --static && $(MAKE) && cp libz.a $(LIB_DIR) && cp zlib.h $(INC_DIR)

$(LIB_DIR)/libbz2.a:
	$(CXX_ENV)cd $(BZIP2_DIR) && $(MAKE) && cp libbz2.a $(LIB_DIR) && cp bzlib.h $(INC_DIR)

$(LIB_DIR)/liblzma.a:
	$(CXX_ENV)cd $(XZ_DIR) && ./configure --disable-xz --disable-xzdec --disable-lzmadec --disable-lzmainfo --disable-lzma-links --disable-scripts --disable-doc --prefix=$(CWD) && make && make install

$(LIB_DIR)/liblz4.a:
	$(CXX_ENV)cd $(LZ4_DIR)/lib && make install prefix=$(CWD)

$(LIB_DIR)/libsnappy.a:
	$(CXX_ENV)cd $(SNAPPY_DIR) && ./configure --prefix=$(CWD) && make && make install

$(LIB_DIR)/libzstd.a:
	$(CXX_ENV)cd $(ZSTD_DIR)/lib && make install prefix=$(CWD)

$(LIB_DIR)/librocksdb.a: $(LIB_DIR)/liblz4.a $(LIB_DIR)/libsnappy.a $(LIB_DIR)/libzstd.a
	$(CXX_ENV)cd $(ROCKSDB_DIR) && make install-static INSTALL_PATH=$(CWD)

$(LIB_DIR)/libpfor.a:
	$(CXX_ENV)cd $(PFOR_DIR) && make && ar rvs libpfor.a *.o && cp *.h $(INC_DIR) && cp libpfor.a $(LIB_DIR)

$(INC_DIR)/tsl/*.h:
	cd $(ROBINMAP_DIR) && cp -r include/tsl $(INC_DIR)

copytsl: $(INC_DIR)/tsl/*.h

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(LIB_DEPS)
	$(CXX_ENV)$(CXX) $(CXX_FLAGS) $(LD_INCLUDE_FLAGS) -c $< -o $@

.prebuild:
	@if [ ! -d $(BIN_DIR) ]; then mkdir -p $(BIN_DIR); fi
	@if [ ! -d $(OBJ_DIR) ]; then mkdir -p $(OBJ_DIR); fi
	@if [ ! -d $(LIB_DIR) ]; then mkdir -p $(LIB_DIR); fi
	@if [ ! -d $(INC_DIR) ]; then mkdir -p $(INC_DIR); fi
	@if [ ! -d $(XZ_DIR) ]; then cd $(DEP_DIR) && tar xzvf xz-5.2.4.tar.gz; fi
	@if [ ! -d $(ZLIB_DIR) ]; then cd $(DEP_DIR) && tar xzvf zlib-1.2.11.tar.gz; fi
	@if [ ! -d $(BZIP2_DIR) ]; then cd $(DEP_DIR) && tar xzvf bzip2-1.0.6.tar.gz; fi
	@if [ ! -d $(LZ4_DIR) ]; then cd $(DEP_DIR) && tar xzvf lz4-1.8.0.tar.gz; fi
	@if [ ! -d $(SNAPPY_DIR) ]; then cd $(DEP_DIR) && tar xzvf snappy-1.1.4.tar.gz; fi
	@if [ ! -d $(ZSTD_DIR) ]; then cd $(DEP_DIR) && tar xzvf zstd-1.4.0.tar.gz; fi
	@if [ ! -d $(ROCKSDB_DIR) ]; then cd $(DEP_DIR) && unzip rocksdb-6.1.2.zip; fi
	@if [ ! -d $(PFOR_DIR) ]; then cd $(DEP_DIR) && unzip TurboPFor-master.zip; fi
	@if [ ! -d $(ROBINMAP_DIR) ]; then cd $(DEP_DIR) && unzip robin-map-0.6.1.zip; fi

# run .pre-build before we make anything at all.
-include .prebuild

clean:
	-$(RM) -f $(OBJ_DIR)/*.o
	-$(RM) -f $(BIN_DIR)/$(EXE)

clean-all:
	$(RM) -r $(BIN_DIR)
	$(RM) -r $(LIB_DIR)
	$(RM) -r $(OBJ_DIR)
	$(RM) -r $(INC_DIR)
	$(RM) -r $(CWD)/share
	$(RM) -r $(CWD)/usr
	cd $(XZ_DIR) && $(MAKE) clean
	cd $(ZLIB_DIR) && $(MAKE) clean
	cd $(BZIP2_DIR) && $(MAKE) clean
	cd $(SNAPPY_DIR) && $(MAKE) clean
	cd $(LZ4_DIR)/lib && $(MAKE) clean
	cd $(ZSTD_DIR)/lib && $(MAKE) clean
	cd $(PFOR_DIR) && $(MAKE) clean && $(RM) libfse.a
	cd $(ROCKSDB_DIR) && $(MAKE) clean
