#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
import warnings

from proteindf_tools import constants
from proteindf_tools.functions import deprecated


class TestFunctions(unittest.TestCase):
    def test_deprecated_function_without_reason(self):
        @deprecated
        def sample_func(a, b):
            """Sample function docstring."""
            return a + b

        # functools.wraps preserves metadata
        self.assertEqual(sample_func.__name__, "sample_func")
        self.assertEqual(sample_func.__doc__, "Sample function docstring.")

        # DeprecationWarning emitted on call
        with warnings.catch_warnings(record=True) as recorded_warnings:
            warnings.simplefilter("always")
            result = sample_func(2, 3)

        self.assertEqual(result, 5)
        self.assertEqual(len(recorded_warnings), 1)
        self.assertTrue(issubclass(recorded_warnings[0].category, DeprecationWarning))
        self.assertEqual(
            str(recorded_warnings[0].message),
            "Call to deprecated function sample_func.",
        )

    def test_deprecated_function_with_reason(self):
        reason = "use new_func instead"

        @deprecated(reason)
        def sample_func_with_reason(x):
            """Sample function with reason docstring."""
            return x * 10

        # functools.wraps preserves metadata
        self.assertEqual(sample_func_with_reason.__name__, "sample_func_with_reason")
        self.assertEqual(sample_func_with_reason.__doc__, "Sample function with reason docstring.")

        # DeprecationWarning emitted on call
        with warnings.catch_warnings(record=True) as recorded_warnings:
            warnings.simplefilter("always")
            result = sample_func_with_reason(4)

        self.assertEqual(result, 40)
        self.assertEqual(len(recorded_warnings), 1)
        self.assertTrue(issubclass(recorded_warnings[0].category, DeprecationWarning))
        self.assertEqual(
            str(recorded_warnings[0].message),
            "Call to deprecated function sample_func_with_reason (use new_func instead).",
        )

    def test_deprecated_class_without_reason(self):
        @deprecated
        class SampleClass:
            """Sample class docstring."""

            def __init__(self, value):
                self.value = value

            def get_value(self):
                return self.value

        # functools.wraps preserves metadata
        self.assertEqual(SampleClass.__name__, "SampleClass")
        self.assertEqual(SampleClass.__doc__, "Sample class docstring.")

        # DeprecationWarning emitted on instantiation
        with warnings.catch_warnings(record=True) as recorded_warnings:
            warnings.simplefilter("always")
            instance = SampleClass(42)

        self.assertEqual(len(recorded_warnings), 1)
        self.assertTrue(issubclass(recorded_warnings[0].category, DeprecationWarning))
        self.assertEqual(
            str(recorded_warnings[0].message),
            "Call to deprecated class SampleClass.",
        )
        self.assertEqual(instance.get_value(), 42)

    def test_deprecated_class_with_reason(self):
        reason = "use NewClass instead"

        @deprecated(reason)
        class SampleClassWithReason:
            """Sample class with reason docstring."""

            def __init__(self, value):
                self.value = value

            def get_value(self):
                return self.value

        # functools.wraps preserves metadata
        self.assertEqual(SampleClassWithReason.__name__, "SampleClassWithReason")
        self.assertEqual(SampleClassWithReason.__doc__, "Sample class with reason docstring.")

        # DeprecationWarning emitted on instantiation
        with warnings.catch_warnings(record=True) as recorded_warnings:
            warnings.simplefilter("always")
            instance = SampleClassWithReason(99)

        self.assertEqual(len(recorded_warnings), 1)
        self.assertTrue(issubclass(recorded_warnings[0].category, DeprecationWarning))
        self.assertEqual(
            str(recorded_warnings[0].message),
            "Call to deprecated class SampleClassWithReason (use NewClass instead).",
        )
        self.assertEqual(instance.get_value(), 99)

    def test_deprecated_invalid_argument(self):
        with self.assertRaises(TypeError):
            deprecated(12345)


class TestConstants(unittest.TestCase):
    def test_au_ang_conversion(self):
        # constants: AU2ANG * ANG2AU == 1
        self.assertEqual(constants.AU2ANG * constants.ANG2AU, 1.0)
        self.assertAlmostEqual(constants.AU2ANG * constants.ANG2AU, 1.0)


if __name__ == "__main__":
    unittest.main()
