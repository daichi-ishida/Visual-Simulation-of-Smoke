#version 330 core

in vec3 texCoords;

uniform sampler3D densityTex;
uniform vec3 eyePos;
uniform vec3 lightPos;
uniform vec3 lightIntensity;
uniform float absorption;
 
out vec4 color;

void main()
{
    // diagonal of the cube
    const float maxDist = sqrt(3.0);
 
    const int numSamples = 128;
    const float scale = maxDist/float(numSamples);
 
    const int numLightSamples = 32;
    const float lscale = maxDist / float(numLightSamples);
 
    // assume all coordinates are in texture space
    vec3 pos = texCoords.xyz;
    vec3 eyeDir = normalize(pos-eyePos)*scale;
 
    // transmittance
    float T = 1.0;
    // in-scattered radiance
    vec3 Lo = vec3(0.0);
 
    for (int i=0; i < numSamples; ++i)
    {
        // sample density
        float density = texture(densityTex, pos).x;
        // tester = vec3(density);
        // skip empty space
        if (density > 0.0)
        {
            // attenuate ray-throughput
            T *= 1.0-density*scale*absorption;
            if (T <= 0.01)
                break;
 
            // point light dir in texture space
            vec3 lightDir = normalize(lightPos-pos)*lscale;
 
            // sample light
            float Tl = 1.0; // transmittance along light ray
            vec3 lpos = pos + lightDir;
 
            for (int s=0; s < numLightSamples; ++s)
            {
                float ld = texture(densityTex, lpos).x;
                Tl *= 1.0-absorption*lscale*ld;
 
                if (Tl <= 0.01)
                    break;
 
                lpos += lightDir;
            }
 
            vec3 Li = lightIntensity*Tl;
 
            Lo += Li*T*density*scale;
        }
 
        pos += eyeDir;
    }
 
    color.xyz = Lo;
    color.w = 1.0-T;
}