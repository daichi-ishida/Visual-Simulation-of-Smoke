#version 330 core

in vec3 texPos;

uniform sampler3D densityTex;
uniform vec3 eyePos;
uniform vec3 lightPos;
uniform vec3 lightIntensity;
uniform float scale;
uniform float absorption;
uniform vec3 ratio;
 
out vec4 color;

void main()
{
    float numSamples = 2.0/scale;
    // assume all coordinates are in texture space
    vec3 pos = 2.0 * texPos - vec3(1.0);
    pos *= ratio;
    vec3 eyeDir = normalize(pos-eyePos)*scale;
 
    // // transmittance
    // float T = 1.0;
    // // in-scattered radiance
    // vec3 Lo = vec3(0.0);
 
    // for (int i=0; i < numSamples; ++i)
    // {
    //     // sample density
    //     float density = texture(densityTex, texPos).x;
    //     // skip empty space
    //     if (density > 0.0)
    //     {
    //         // attenuate ray-throughput
    //         T *= 1.0 - density*scale*absorption;
    //         if (T <= 0.01)
    //         {
    //             break;
    //         }
    //         // point light dir in texture space
    //         vec3 lightDir = normalize(lightPos-pos)*scale;
 
    //         // sample light
    //         float Tl = 1.0; // transmittance along light ray
    //         vec3 lpos = texPos + lightDir;
 
    //         for (int s=0; s < numSamples; ++s)
    //         {
    //             float ld = texture(densityTex, lpos).x;
    //             Tl *= 1.0-absorption*scale*ld;
 
    //             if (Tl <= 0.01)
    //                 break;
 
    //             lpos += lightDir;
    //         }
 
    //         vec3 Li = lightIntensity*Tl;
 
    //         Lo += Li*T*density*scale;
    //     }
 
    //     pos += eyeDir;
    // }
 
    color.xyz = pos;
    color.w = 1.0 - T;
}