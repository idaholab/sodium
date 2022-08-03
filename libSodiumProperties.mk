#
# build libSodiumProperties using libMesh's build system
#
LIBSODIUM_DIR       := $(SODIUM_DIR)/contrib/libSodiumProperties

LIBSODIUM_srcfiles  :=
LIBSODIUM_srcfiles  += $(LIBSODIUM_DIR)/Na_Golden.cpp

LIBSODIUM_objects   := $(patsubst %.cpp, %.$(obj-suffix), $(LIBSODIUM_srcfiles))
LIBSODIUM_deps      := $(patsubst %.$(obj-suffix), %.$(obj-suffix).d, $(LIBSODIUM_objects))
LIBSODIUM_LIB       := $(LIBSODIUM_DIR)/libSodiumProperties-$(METHOD).la

app_INCLUDES += -I$(SODIUM_DIR)
app_LIBS += $(LIBSODIUM_LIB)

$(LIBSODIUM_LIB): $(LIBSODIUM_objects)
	@echo "Linking Library "$@"..."
	@$(libmesh_LIBTOOL) --tag=CC $(LIBTOOLFLAGS) --mode=link --quiet \
	  $(libmesh_CC) $(libmesh_CFLAGS) -o $@ $(LIBSODIUM_objects) $(libmesh_LDFLAGS) $(EXTERNAL_FLAGS) -rpath $(LIBSODIUM_DIR)
	@$(libmesh_LIBTOOL) --mode=install --quiet install -c $(LIBSODIUM_LIB) $(LIBSODIUM_DIR)

$(app_EXEC): $(LIBSODIUM_LIB)

-include $(LIBSODIUM_deps)

cleanlibsodium:
	@echo "Cleaning libSodiumProperties"
	@rm -f $(LIBSODIUM_objects)
	@rm -f $(LIBSODIUM_deps)
	@rm -f $(LIBSODIUM_LIB)
	@rm -f $(LIBSODIUM_DIR)/libSodiumProperties-$(METHOD)*.dylib
	@rm -f $(LIBSODIUM_DIR)/libSodiumProperties-$(METHOD)*.so*
	@rm -f $(LIBSODIUM_DIR)/libSodiumProperties-$(METHOD)*.a
