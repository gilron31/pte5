// The one macro that lets a value type be shared by the host and the device.
//
// ---------------------------------------------------------------------------
// Why this exists, and why it does not break the extension rule
// ---------------------------------------------------------------------------
//
// gpu/core/types.h states the rule this project builds on:
//
//   .cu/.cuh iff the file contains device code, because one unguarded
//   __device__ function in a header promotes every translation unit that
//   includes it, transitively, to nvcc.
//
// `PTE_HOST_DEVICE` is the sanctioned escape hatch, and it keeps that rule
// intact rather than bending it: expanded by g++ it is *empty*, so a host
// translation unit including core/u256.h sees ordinary inline functions and
// never learns that CUDA exists. Expanded by nvcc it marks the same functions
// as callable from a kernel. Nothing is promoted to nvcc by inclusion.
//
// This is the same bargain device_append_buffer.h already makes with its
// `#ifdef __CUDACC__` guard around `push`, and it carries the same caveat,
// restated here because this header multiplies it across a whole type:
//
//   **the two definitions must be identical apart from the annotation.**
//
// Formally the guard gives a function a different declaration in host and
// device translation units, which is an ODR violation. It is benign in exactly
// the way the existing one is -- the bodies are the same source text, the types
// are trivially copyable aggregates whose layout the macro cannot touch, and
// the device copy is emitted only into device code. It stops being benign the
// moment someone writes `#ifdef __CUDA_ARCH__` around something that *changes a
// value*. Dispatching on __CUDA_ARCH__ to reach a faster instruction for the
// same arithmetic result (u256.h's 64x64 multiply does exactly this) is fine
// and is tested for agreement in gpu/bigint/u256_device_test.cu; branching to a
// different answer is not.
#pragma once

#if defined(__CUDACC__)
#define PTE_HOST_DEVICE __host__ __device__
#else
#define PTE_HOST_DEVICE
#endif
