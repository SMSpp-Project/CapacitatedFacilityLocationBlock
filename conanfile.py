from conans import ConanFile, CMake, tools


class SmsppConan(ConanFile):
    name = "mcfblock"
    version = "0.4.0"
    description = "SMS++ Block and Solver for linear Min-Cost Flow problems"
    topics = ("conan", "smspp", "mcfblock")
    url = "https://gitlab.com/smspp/mcfblock"
    homepage = "https://gitlab.com/smspp/mcfblock"
    license = "GPL-3.0-only"
    generators = "cmake"

    settings = "os", "arch", "compiler", "build_type"
    options = {"shared": [True, False], "fPIC": [True, False]}
    default_options = {"shared": False, "fPIC": True}

    requires = (
        "mcfclass/20200227@smspp/testing",
        "smspp/0.4.0@smspp/testing"
    )

    exports_sources = [
        "CMakeLists.txt",
        "src/*",
        "include/*",
        "cmake/*",
        "tools/*"
    ]

    def source(self):
        tools.replace_in_file(
            "CMakeLists.txt",
            '''LANGUAGES C CXX)''',
            '''LANGUAGES C CXX)\n''' +
            '''include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)\n''' +
            '''conan_basic_setup()'''
        )

    def _configure_cmake(self):
        cmake = CMake(self)
        cmake.definitions["BUILD_TESTING"] = False
        cmake.configure()
        return cmake

    def build(self):
        cmake = self._configure_cmake()
        cmake.build()

    def package(self):
        self.copy(pattern="LICENSE", dst="licenses")
        cmake = self._configure_cmake()
        cmake.install()

    def package_info(self):
        self.cpp_info.includedirs = ["include", "include/SMS++"]
        self.cpp_info.libs = ["MCFBlock"]
