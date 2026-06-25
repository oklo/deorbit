// Toolchain verification: create the GPU device, load hello.metallib, dispatch
// vadd over unified-memory buffers, read back. If c[]==3, the whole Metal path
// (compiler + metal-cpp host + dispatch + shared memory) works.
#define NS_PRIVATE_IMPLEMENTATION
#define MTL_PRIVATE_IMPLEMENTATION
#define CA_PRIVATE_IMPLEMENTATION
#include <Metal/Metal.hpp>
#include <Foundation/Foundation.hpp>
#include <cstdio>
#include <vector>

int main() {
    MTL::Device* dev = MTL::CreateSystemDefaultDevice();
    if (!dev) { printf("no Metal device\n"); return 1; }
    printf("GPU: %s\n", dev->name()->utf8String());

    NS::Error* err = nullptr;
    auto libpath = NS::String::string("hello.metallib", NS::UTF8StringEncoding);
    auto lib = dev->newLibrary(libpath, &err);
    if (!lib) { printf("library load failed: %s\n",
                       err ? err->localizedDescription()->utf8String() : "?"); return 1; }
    auto fn = lib->newFunction(NS::String::string("vadd", NS::UTF8StringEncoding));
    auto pso = dev->newComputePipelineState(fn, &err);
    auto queue = dev->newCommandQueue();

    const int n = 1 << 20;                       // 1M elements
    std::vector<float> a(n, 1.0f), b(n, 2.0f);
    auto ba = dev->newBuffer(a.data(), n * sizeof(float), MTL::ResourceStorageModeShared);
    auto bb = dev->newBuffer(b.data(), n * sizeof(float), MTL::ResourceStorageModeShared);
    auto bc = dev->newBuffer(n * sizeof(float), MTL::ResourceStorageModeShared);

    auto cmd = queue->commandBuffer();
    auto enc = cmd->computeCommandEncoder();
    enc->setComputePipelineState(pso);
    enc->setBuffer(ba, 0, 0); enc->setBuffer(bb, 0, 1); enc->setBuffer(bc, 0, 2);
    enc->dispatchThreads(MTL::Size(n, 1, 1), MTL::Size(256, 1, 1));
    enc->endEncoding(); cmd->commit(); cmd->waitUntilCompleted();

    float* c = (float*)bc->contents();
    printf("c[0]=%.1f c[%d]=%.1f  (expect 3.0)  -> %s\n",
           c[0], n - 1, c[n - 1], (c[0] == 3.0f && c[n - 1] == 3.0f) ? "OK" : "FAIL");
    return 0;
}
