// Scene shaders (Unit 25a, VR5, VR8): triangles, opaque or into the weighted blended order-independent
// transparency targets, and lines as screen-space quads. Compiled at build time with dxc (VQ2).

cbuffer Frame : register(b0) {
    row_major float4x4 view_projection;
    float3 light;        // unit, toward the light, world space
    float ambient;
    float2 viewport;     // target size in pixels
    uint flags;          // bit 0 lit, bit 1 field, bit 2 transparent pass, bits 3-4 finish (Unit 26g)
    float line_width;    // pixels at the target's resolution
    float4 colormap[256];
    float4 no_data;
    // The enhanced look (Unit 26f, KD2, KD4); zero for classic.
    float3 eye;          // unit, toward the viewer (orthographic), or the camera position (perspective)
    float style;         // 0 classic, 1 enhanced
    float3 up;           // unit, the camera's up, world space
    float perspective;   // 1 when eye is a position
    float3 fill;         // unit, toward the fill light, world space
    float spare;
    float4 bg_top;
    float4 bg_bottom;
    float4 view_rect;    // left, top, width, height of the view in target pixels
    float3 half_vector;  // orthographic: normalize(light + eye), the same at every pixel
    float spare2;
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
    float3 world : WORLD;  // um, for the enhanced look's eye direction in perspective
};

MeshOut mesh_vs(MeshIn v) {
    MeshOut o;
    o.position = mul(float4(v.position, 1.0), view_projection);
    o.normal = v.normal;
    o.value = v.value;
    o.color = v.color;
    o.world = v.position;
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
    if ((flags & 1u) && style > 0.5) {
        // Enhanced: key, fill and rim lights, a sky/ground hemisphere ambient and a soft highlight
        // (reduced on field maps, KQ6); two-sided by turning the normal toward the eye.
        float3 n = normalize(p.normal);
        const float3 v = perspective > 0.5 ? normalize(eye - p.world) : eye;
        if (dot(n, v) < 0.0) n = -n;
        const uint finish = (flags >> 3) & 3u;  // Unit 26g (MD1): 0 matte, 1 metal, 2 glass
        const bool metal = finish == 1u;
        const float key = saturate(dot(n, light));
        const float fil = saturate(dot(n, fill));
        // Metal (MD2): the sky/ground ambient along the reflected eye ray, as a fake environment, so faces
        // turned up read bright and faces turned down read dark; aluminium as its base colour.
        const float3 env = metal ? reflect(-v, n) : n;
        const float hemi = lerp(0.30, 0.46, 0.5 + 0.5 * dot(env, up));
        const float lum = hemi + 0.60 * key + 0.18 * fil;
        const float field = (flags & 2u) ? 1.0 : 0.0;
        if (metal) c.rgb = float3(0.70, 0.72, 0.76);
        // Powers by multiplication (s^40 = s^32 s^8, s^64, r^3), no pow: WARP runs this per pixel on the CPU (KG6).
        const float3 h = perspective > 0.5 ? normalize(light + v) : half_vector;
        const float s = saturate(dot(n, h));
        const float s2 = s * s, s4 = s2 * s2, s8 = s4 * s4, s32 = s8 * s8 * s8 * s8;
        const float r = 1.0 - saturate(dot(n, v));
        const float r3 = r * r * r;
        const float spec = metal ? 0.55 * (s32 * s32) : lerp(0.22, 0.08, field) * (s32 * s8);
        const float rim = lerp(0.14, 0.06, field) * r3;
        c.rgb = c.rgb * lum + spec + rim;
        if (finish == 2u) {
            // Glass (MD3): more opaque and lighter toward grazing angles (Schlick-like), as clear as
            // before face on, so the oxide's outline shows where it is seen edge-on.
            const float g = 0.55 * r3;
            c.a = c.a + (1.0 - c.a) * g;
            c.rgb = lerp(c.rgb, float3(1.0, 1.0, 1.0), g);
        }
    } else if (flags & 1u) {
        // Two-sided: a slice or an interface is lit from whichever side the eye is on.
        const float d = abs(dot(normalize(p.normal), light));
        c.rgb *= ambient + (1.0 - ambient) * d;
    }
    return c;
}

// The enhanced look's background (KD4): a vertical gradient over the view, drawn before the scene.
float4 background_ps(float4 p : SV_Position) : SV_Target {
    const float t = saturate((p.y - view_rect.y) / max(view_rect.w, 1.0));
    return float4(lerp(bg_top.rgb, bg_bottom.rgb, t), 1.0);
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
