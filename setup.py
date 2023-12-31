#!/usr/bin/env python3

import os
import sys
import glob
import shlex
import shutil
from setuptools import setup, find_packages, Extension, Command
from setuptools.command.test import test
from setuptools.command.build_ext import build_ext
from setuptools.command.install import install
import subprocess

setup_src = os.path.dirname(os.path.realpath(__file__))


class CMakeExtension(Extension):
    """Initialise the name of a CMake extension."""

    def __init__(self, name):
        super().__init__(name, sources=[])


class CMakeBuild(build_ext):
    """Build and configure a CMake extension."""

    def run(self):
        link_args = []
        if sys.platform == "darwin":
            link_args.append("-Wl,-rpath,@loader_path")
        for ext in self.extensions:
            ext.link_args = link_args
            self.build_cmake(ext)
        super().run()

    def build_cmake(self, ext):
        src = os.path.join(setup_src, "pyhub", "libs")

        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        # Ajout de la recherche automatique du compilateur Fortran
        try:
            fortran_compiler = os.environ["FC"]
        except KeyError:
            # Si FC n'est pas défini, essayez de trouver le compilateur Fortran par défaut
            try:
                fortran_compiler = (
                    os.popen("which gfortran").read().strip() or
                    os.popen("which ifort").read().strip()
                )
            except FileNotFoundError:
                raise EnvironmentError("Compilateur Fortran non trouvé. Définissez la variable d'environnement FC ou assurez-vous que gfortran ou ifort est installé.")

        cmake_args = [f"-S{src}", f"-B{self.build_temp}", f"-DCMAKE_Fortran_COMPILER={fortran_compiler}"]
#        cmake_args = [f"-S{src}", f"-B{self.build_temp}"]
        if os.getenv("CMAKE_CONFIGURE_ARGS"):
            cmake_args += os.getenv("CMAKE_CONFIGURE_ARGS").split()

        self.announce("Configuring")
        self.spawn(["cmake", *cmake_args])

        build_args = []
        if os.getenv("CMAKE_BUILD_ARGS"):
            cmake_args += os.getenv("CMAKE_BUILD_ARGS").split()
        if getattr(self, "parallel", False):
            build_args.append(f"-j{self.parallel}")

        self.announce("Building")
        self.spawn(["cmake", "--build", self.build_temp, *build_args])


class CleanCommand(Command):
    """Clean up files resulting from compilation except for .so shared objects."""

    CLEAN_FILES = ["build", "dist", "*.egg-info"]
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        for path_spec in self.CLEAN_FILES:
            paths = glob.glob(os.path.normpath(os.path.join(setup_src, path_spec)))
            for path in paths:
                if not str(path).startswith(setup_src):
                    # In case CLEAN_FILES contains an address outside the package
                    raise ValueError("%s is not a path inside %s" % (path, setup_src))
                shutil.rmtree(path)


class DiscoverTests(test):
    """Discover and dispatch tests."""

    user_options = [
        ("include-veryslow", "v", "Include tests marked as veryslow"),
        ("include-slow", "s", "Include tests marked as slow"),
        ("pytest-args=", "p", "Extra arguments for pytest"),
    ]

    def initialize_options(self):
        test.initialize_options(self)
        self.include_veryslow = False
        self.include_slow = True
        self.pytest_args = ""

    def finalize_options(self):
        pass

    def run_tests(self):
        # Only import pytest in this scope
        import pytest

        src = os.path.join(setup_src, "pyhub", "tests")

        test_args = []
        if not (self.include_slow and self.include_veryslow):
            test_args.append("-m not (slow or veryslow)")
        elif not self.include_veryslow:
            test_args.append("-m not veryslow")
        elif not self.include_slow:
            test_args.append("-m not slow")
        test_args += shlex.split(self.pytest_args)

        pytest.main([src, *test_args])

class InstallCommand(install):
    def run(self):
        # Exécute la méthode d'installation de la classe parente
        install.run(self)

        # Installez vos dépendances ici si elles ne sont pas déjà présentes
        dependencies = [
            "numpy>=1.19.0",
            "scipy>=1.1.0",
            "h5py>=2.7",
        ]
        for dependency in dependencies:
            try:
                __import__(dependency)
            except ImportError:
                subprocess.call(['pip', 'install', dependency])

from distutils.command.build import build

build.sub_commands = [c for c in build.sub_commands if c[0] == "build_ext"] + [
    c for c in build.sub_commands if c[0] != "build_ext"
]

setup(
    packages=find_packages(exclude=["*examples*"]),
    include_package_data=True,
    install_requires=[
    "numpy>=1.19.0",
    "scipy>=1.1.0",
    "h5py>=2.7",
    ],
    ext_modules=[CMakeExtension("pyhub/libs")],
    cmdclass={
        "build_ext": CMakeBuild,
        "test": DiscoverTests,
        "clean": CleanCommand,
        'install': InstallCommand,
    },
    zip_safe=False,
)