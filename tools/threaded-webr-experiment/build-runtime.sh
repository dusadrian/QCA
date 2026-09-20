#!/usr/bin/env bash
set -euo pipefail

experiment_dir="$(cd "$(dirname "$0")" && pwd)"
qca_dir="$(cd "$experiment_dir/../.." && pwd)"
build_dir="$experiment_dir/build"
emscripten_cache="$build_dir/emscripten-cache"
image="${WEBR_IMAGE:-ghcr.io/r-wasm/webr:v0.6.0}"

mkdir -p "$build_dir" "$emscripten_cache"
rm -rf "$build_dir/runtime" "$build_dir/repo"
rm -f "$build_dir"/QCA_*.tar.gz "$build_dir"/QCA_*.tgz

(cd "$build_dir" && R CMD build --no-build-vignettes "$qca_dir")
qca_tarball="$(find "$build_dir" -maxdepth 1 -name 'QCA_*.tar.gz' -print -quit)"
if [[ -z "$qca_tarball" ]]; then
    echo "QCA source tarball was not created." >&2
    exit 1
fi

docker image inspect "$image" --format '{{.Id}}' > "$build_dir/docker-image-id"

if [[ ! -f "$emscripten_cache/.qca-seeded" ]]; then
    docker run --rm --platform linux/amd64 \
        -v "$emscripten_cache:/cache" \
        "$image" \
        sh -lc 'cp -a /opt/emsdk/upstream/emscripten/cache/. /cache/ && touch /cache/.qca-seeded'
fi

docker run --rm --platform linux/amd64 \
    -v "$experiment_dir:/experiment" \
    -v "$emscripten_cache:/emscripten-cache" \
    -e EM_CACHE=/emscripten-cache \
    -w /experiment \
    "$image" \
    sh -lc '
        set -eu
        export EMCC_CORES="${WEBR_BUILD_JOBS:-4}"
        patch -d /opt/emsdk/upstream/emscripten -p1 \
            < /experiment/emscripten-settings-file.patch
        rm -rf /tmp/threaded-webr
        cp -a /opt/webr /tmp/threaded-webr
        cd /tmp/threaded-webr
        patch -p1 < /experiment/webr-pthreads.patch

        trap '\''
            mkdir -p /experiment/build/diagnostics
            find R/build -path "*/build/config.log" \
                -exec cp {} /experiment/build/diagnostics/config.log \; 2>/dev/null || true
        '\'' EXIT

        version="$(cat R/R-VERSION)"
        r_build="R/build/R-${version}/build"

        grep -Ilr /opt/webr "$r_build" | while IFS= read -r file; do
            sed -i "s#/opt/webr#/tmp/threaded-webr#g" "$file"
        done
        for makeconf in \
            "$r_build/Makeconf" \
            "$r_build/etc/Makeconf" \
            "wasm/R-${version}/lib/R/etc/Makeconf"; do
            sed -i \
                -e "s#^CPPFLAGS =#CPPFLAGS = -pthread#" \
                -e "s#^CFLAGS =#CFLAGS = -pthread#" \
                -e "s#^LDFLAGS =#LDFLAGS = -pthread#" \
                -e "s#^SHLIB_LDFLAGS =#SHLIB_LDFLAGS = -pthread#" \
                "$makeconf"
        done

        make -C "$r_build" clean
        sed -i "s/shlib  cairodevice/shlib/" \
            "$r_build/src/library/grDevices/src/Makefile"
        rm -f "R/build/state/R-${version}/r-stage2"

        # The stock image Cairo dependencies were compiled without shared
        # memory.  Skip only that optional backend; rebuild the core and the rest
        # of the WebR base package payload coherently with pthread support.
        make NPROCS="${WEBR_BUILD_JOBS:-4}" webr

        mkdir -p /experiment/build/runtime
        cp -a dist/. /experiment/build/runtime/
        cp "$r_build/src/main/R.wasm" /experiment/build/runtime/R.wasm
        cp "$r_build/lib/libRblas.so" /experiment/build/runtime/libRblas.so
        cp "$r_build/lib/libRlapack.so" /experiment/build/runtime/libRlapack.so
        cat "wasm/R-${version}/pre.js" "$r_build/src/main/R.bin" \
            > /experiment/build/runtime/R.js
        node /experiment/patch-r-runtime.mjs \
            /experiment/build/runtime/R.js
        node /experiment/patch-webr-worker.mjs \
            /experiment/build/runtime/webr-worker.js

        rm -rf /tmp/qca-bin
        qca_tarball="$(find /experiment/build -maxdepth 1 -name "QCA_*.tar.gz" -print -quit)"
        mkdir -p /tmp/qca-bin
        Rscript -e "
            suppressPackageStartupMessages(library(rwasm))
            options(rwasm.webr_root = \"/tmp/threaded-webr\")
            rwasm:::wasm_build(
                \"QCA\",
                \"${qca_tarball}\",
                \"/tmp/qca-bin\",
                TRUE
            )
        "
        mkdir -p /experiment/build/repo/bin/emscripten/contrib/4.6
        cp -a /tmp/qca-bin/. \
            /experiment/build/repo/bin/emscripten/contrib/4.6/
    '

qca_binary="$(find "$build_dir/repo/bin/emscripten/contrib" -name 'QCA_*.tgz' -print -quit)"
if [[ -z "$qca_binary" ]]; then
    echo "Threaded QCA WebAssembly package was not created." >&2
    exit 1
fi

mkdir -p "$build_dir/qca-package"
tar -xzf "$qca_binary" -C "$build_dir/qca-package"
qca_module="$(find "$build_dir/qca-package" -path '*/libs/QCA.so' -print -quit)"
if [[ -z "$qca_module" ]]; then
    echo "QCA WebAssembly side module was not found in the binary package." >&2
    exit 1
fi
cp "$qca_module" "$build_dir/QCA.so"

printf '%s\n' "Threaded WebR runtime: $build_dir/runtime"
printf '%s\n' "Threaded QCA side module: $build_dir/QCA.so"
