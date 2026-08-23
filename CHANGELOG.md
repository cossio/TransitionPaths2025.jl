# Changelog

All notable changes to this project will be documented in this file. The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

- Declare the test environment as a Pkg `[workspace]` member, so test deps resolve in a single shared manifest.
- Require Julia 1.12 (needed for `[workspace]` support); CI now tests on 1.12.

## v0.1.0

- Init package.