# TOMBO_TEXT_OUTPUT Module

## Container Dependencies

This module requires both `ont-tombo` and `pandas` to function properly. The module uses Wave (Seqera Containers) for automatic multi-package container building following nf-core best practices.

### Recommended Usage with Wave (Modern nf-core Standard)

**Use the Wave profile** (Automatic container building):
```bash
nextflow run main.nf -profile wave
```

Wave will automatically:
- Build a container from the `environment.yml` specification
- Include both ont-tombo and pandas
- Cache the container for future runs
- Support multiple architectures (AMD64, ARM64)

### Alternative Usage Options:

1. **Conda Profile** (Local installation):
   ```bash
   nextflow run main.nf -profile conda
   ```
   Installs packages locally from `environment.yml`.

2. **Enable Wave in config** (For custom setups):
   Add to your `nextflow.config` or command line:
   ```bash
   nextflow run main.nf -profile docker -with-wave
   ```

### Environment Specification

The `environment.yml` follows nf-core standards:
- Includes Renovate comments for automated dependency updates
- Specifies exact versions for reproducibility
- Uses conda-forge and bioconda channels

### Testing

The nf-test suite includes:
- **Stub mode test**: ✓ Always passes (validates module structure)
- **Full execution test**: ✓ Passes with `-profile wave` or `-profile conda`

### Wave Benefits

Wave provides:
- **Automatic container building** from conda specifications
- **Multi-architecture support** (AMD64, ARM64, etc.)
- **Container caching** for faster subsequent runs
- **No manual container building required**
- **Automatic dependency resolution**
