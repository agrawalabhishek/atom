# Copyright (c) 2014-2016 Kartik Kumar, Dinamica Srl (me@kartikkumar.com)
# Distributed under the MIT License.
# See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT

# Set project main file.
set(MAIN_SRC
  "${SRC_PATH}/main.cpp"
)

# Set project source files.
#set(SRC
# "${SRC_PATH}/tools.cpp"
#)

# Set project test source files.
set(TEST_SRC
  "${TEST_SRC_PATH}/testAtom.cpp"
#  "${TEST_SRC_PATH}/testAtomSolver.cpp"
  "${TEST_SRC_PATH}/testConvertCartesianStateToTwoLineElements.cpp"
  "${TEST_SRC_PATH}/testPrintFunctions.cpp"
)
