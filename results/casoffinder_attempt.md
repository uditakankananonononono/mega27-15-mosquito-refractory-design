# Cas-OFFinder local run attempt (blocked, documented)

Cas-OFFinder 2.4.1 linux x86-64 binary (GitHub release asset cas-offinder_linux_x86-64.zip)
fails at startup in this sandbox: `clGetPlatformIDs Failed: -1001` - the binary requires an
OpenCL platform (GPU or CPU-emulated such as POCL), and none is available in this environment.

Consequence: the Appendix M throughput comparison remains against the PUBLISHED Cas-OFFinder
numbers (Bae, Park, Kim 2014), and is labeled as such there. A first-party local run stays
queued behind an OpenCL runtime. This file exists so the attempt is on record rather than
silently dropped.
