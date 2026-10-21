// Full-screen passes (Unit 25a): composite the order-independent transparency targets over the
// opaque image, and average the supersampled image 2 x 2 into the output.

Texture2D<float4> accum_texture : register(t0);
Texture2D<float> reveal_texture : register(t1);
Texture2D<float4> source_texture : register(t0);

float4 fullscreen_vs(uint id : SV_VertexID) : SV_Position {
    const float2 uv = float2((id << 1) & 2, id & 2);
    return float4(uv * float2(2.0, -2.0) + float2(-1.0, 1.0), 0.0, 1.0);
}

// Blended over the opaque image with (src alpha, 1 - src alpha).
float4 composite_ps(float4 p : SV_Position) : SV_Target {
    const int2 xy = int2(p.xy);
    const float reveal = reveal_texture.Load(int3(xy, 0));
    if (reveal >= 1.0) discard;
    const float4 a = accum_texture.Load(int3(xy, 0));
    return float4(a.rgb / max(a.a, 1e-5), 1.0 - reveal);
}

float4 resolve_ps(float4 p : SV_Position) : SV_Target {
    const int2 xy = int2(p.xy) * 2;
    return 0.25 * (source_texture.Load(int3(xy, 0)) + source_texture.Load(int3(xy + int2(1, 0), 0)) +
                   source_texture.Load(int3(xy + int2(0, 1), 0)) + source_texture.Load(int3(xy + int2(1, 1), 0)));
}
