// Scene shaders (Unit 25a, VR5, VR8): triangles, opaque or into the weighted blended order-independent
// transparency targets, and lines as screen-space quads. Compiled at build time with dxc (VQ2).

cbuffer Frame : register(b0) {
    row_major float4x4 view_projection;
    float3 light;        // unit, toward the light, world space
    float ambient;
    float2 viewport;     // target size in pixels
    uint flags;          // bit 0 lit, bit 1 field, bit 2 transparent pass
    float line_width;    // pixels at the target's resolution
    float4 colormap[256];
    float4 no_data;
};

struct MeshIn {
    float3 position : POSITION;
    float3 normal : NORMAL;
    float value : VALUE;
    float4 color : COLOR;
};

struct MeshOut {
    float4 position : SV_Position;
    float3 normal : NORMAL;
    float value : VALUE;
    float4 color : COLOR;
};

MeshOut mesh_vs(MeshIn v) {
    MeshOut o;
    o.position = mul(float4(v.position, 1.0), view_projection);
    o.normal = v.normal;
    o.value = v.value;
    o.color = v.color;
    return o;
}

float4 base_color(MeshOut p) {
    float4 c = p.color;
    if (flags & 2u) {
        if (isnan(p.value)) {
            c = no_data;
        } else {
            const float x = saturate(p.value) * 255.0;
            const uint i = min((uint)x, 254u);
            c = float4(lerp(colormap[i].rgb, colormap[i + 1].rgb, x - (float)i), p.color.a);
        }
    }
    if (flags & 1u) {
        // Two-sided: a slice or an interface is lit from whichever side the eye is on.
        const float d = abs(dot(normalize(p.normal), light));
        c.rgb *= ambient + (1.0 - ambient) * d;
    }
    return c;
}

float4 opaque_ps(MeshOut p) : SV_Target {
    float4 c = base_color(p);
    c.a = 1.0;
    return c;
}

// Weighted blended OIT: accumulate premultiplied color times a depth weight, and the product of
// (1 - alpha) as revealage.
struct OitOut {
    float4 accum : SV_Target0;
    float reveal : SV_Target1;
};

OitOut transparent_ps(MeshOut p) {
    const float4 c = base_color(p);
    const float z = p.position.z;
    const float w = clamp(pow(min(1.0, c.a * 10.0) + 0.01, 3.0) * 1e8 * pow(1.0 - z * 0.9, 3.0), 1e-2, 3e3);
    OitOut o;
    o.accum = float4(c.rgb * c.a, c.a) * w;
    o.reveal = c.a;
    return o;
}

// Lines: each segment is six vertices carrying both end points; the vertex shader offsets this end
// across the segment by half the width in pixels.
struct LineIn {
    float3 position : POSITION;
    float3 other : OTHER;
    float side : SIDE;      // -1 or +1 across the segment
    float4 color : COLOR;
};

struct LineOut {
    float4 position : SV_Position;
    float4 color : COLOR;
};

LineOut line_vs(LineIn v) {
    const float4 a = mul(float4(v.position, 1.0), view_projection);
    const float4 b = mul(float4(v.other, 1.0), view_projection);
    const float2 sa = a.xy / a.w * viewport * 0.5;
    const float2 sb = b.xy / b.w * viewport * 0.5;
    float2 d = sb - sa;
    const float len = length(d);
    d = len > 0.0 ? d / len : float2(1.0, 0.0);
    const float2 n = float2(-d.y, d.x);
    LineOut o;
    o.position = a;
    o.position.xy += n * v.side * line_width / viewport * a.w;
    o.position.z -= 1e-4 * a.w;  // drawn over the surface it lies on
    o.color = v.color;
    return o;
}

float4 line_ps(LineOut p) : SV_Target {
    return p.color;
}
