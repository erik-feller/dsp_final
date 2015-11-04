from distutils.core import setup, Extension

setup(name="usrp",
      version="0.1",
      ext_modules=[Extension("usrp",["usrp.cc"],
                             include_dirs = ["/usr/local/include"],
                             library_dirs = ["/usr/local/lib"],
                             libraries = ["usrp"],
                             )])

