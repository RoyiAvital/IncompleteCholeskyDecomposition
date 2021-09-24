
set_project("IncompleteCholeskyDecomposition")
set_version("0.0.1")
set_xmakever("2.3.7")

-- set_toolchains("mingw")
add_rules("mode.release")

-- import("core.tool.compiler")



target("IncompleteCholeskyDecomposition")

	-- local flags = compiler.compflags("convolution.cpp")
	-- for _, flag in ipairs(flags) do
		-- print(flag)
	-- end

	set_kind("shared")
	set_languages("c17")
	-- add_files("src/*.cpp", "src/*.c", "src/Image_structures/*.cpp", "src/Patch_match/*.cpp", "src/Reconstruction/*.cpp")
	add_files("IncompleteCholeskyDecomposition.c")
	add_cxflags("-march=native") -- C / C++
	-- add_cxxflags("-g", "-O2", "-DDEBUG") -- C++ Only
	set_optimize("aggressive")
	add_vectorexts("avx2")
	--add_defines("DEBUG")
	--add_linkdirs("$(buildir)/lib")
	--add_links("png", "z")
	-- add_syslinks("png", "z")
	set_warnings("everything")
	
target_end()


-- For configuration run (MinGW):
-- xmake f -p mingw --mingw=D:\Applications\Programming\MinGW -v
-- For compilation
-- xmake -v -w