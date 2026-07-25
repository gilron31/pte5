// pybind11 type caster mapping GMP's mpz_class to and from Python's built-in
// int, so bound functions can take and return native Python integers of any
// size with no precision loss.
//
// The conversion goes through a decimal string. That is not the fastest route,
// but it is exact, needs no knowledge of CPython's internal digit layout, and
// stays correct across CPython versions.
#ifndef PTE_BINDINGS_MPZ_CASTER_H_
#define PTE_BINDINGS_MPZ_CASTER_H_

#include <gmpxx.h>
#include <pybind11/pybind11.h>

#include <string>

namespace pybind11 {
namespace detail {

template <> struct type_caster<mpz_class> {
  PYBIND11_TYPE_CASTER(mpz_class, const_name("int"));

  // Python -> mpz_class. Accepts anything with __index__ (int, numpy integers),
  // and rejects float and str, which is the behaviour a caller expects from an
  // integer parameter.
  bool load(handle src, bool /*convert*/) {
    if (!src)
      return false;
    object index = reinterpret_steal<object>(PyNumber_Index(src.ptr()));
    if (!index) {
      PyErr_Clear();
      return false;
    }
    object text = reinterpret_steal<object>(PyObject_Str(index.ptr()));
    if (!text) {
      PyErr_Clear();
      return false;
    }
    const char *digits = PyUnicode_AsUTF8(text.ptr());
    if (digits == nullptr) {
      PyErr_Clear();
      return false;
    }
    return value.set_str(digits, 10) == 0;
  }

  // mpz_class -> Python int.
  static handle cast(const mpz_class &src, return_value_policy /*policy*/,
                     handle /*parent*/) {
    const std::string digits = src.get_str();
    PyObject *result = PyLong_FromString(digits.c_str(), nullptr, 10);
    if (result == nullptr)
      throw error_already_set();
    return reinterpret_steal<object>(result).release();
  }
};

} // namespace detail
} // namespace pybind11

#endif // PTE_BINDINGS_MPZ_CASTER_H_
