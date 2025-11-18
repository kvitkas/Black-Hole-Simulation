#version 120

uniform vec2 u_resolution;
uniform float u_time;
uniform vec2 u_mouse;

const float PI = 3.14159265359;
const int MAX_STEPS = 256;
const float MAX_DIST = 30.0;
const float STEP_SIZE = 0.04;
const float M = 1.0;

float hash(float n) { return fract(sin(n) * 43758.5453123); }
float noise(vec3 x) {
    vec3 p = floor(x);
    vec3 f = fract(x);
    f = f * f * (3.0 - 2.0 * f);
    float n = p.x + p.y * 57.0 + 113.0 * p.z;
    return mix(mix(mix(hash(n + 0.0), hash(n + 1.0), f.x),
                   mix(hash(n + 57.0), hash(n + 58.0), f.x), f.y),
               mix(mix(hash(n + 113.0), hash(n + 114.0), f.x),
                   mix(hash(n + 170.0), hash(n + 171.0), f.x), f.y), f.z);
}

float fbm(vec3 p) {
    float f = 0.0;
    float amp = 0.5;
    for(int i=0; i<5; i++) {
        f += amp * noise(p);
        p *= 2.02;
        amp *= 0.5;
    }
    return f;
}

vec3 rotateY(vec3 v, float angle) {
    float c = cos(angle);
    float s = sin(angle);
    return vec3(c * v.x + s * v.z, v.y, -s * v.x + c * v.z);
}

vec3 getDiskColor(float r, float phi) {
    if (r < 2.6 || r > 14.0) return vec3(0.0);
    
    float temp = 1.0 / (r - 2.0);
    float n = fbm(vec3(r * 1.5, phi * 4.0, u_time * 0.2));
    
    float intensity = temp * (0.2 + 0.8 * n);
    intensity = pow(intensity, 1.5) * 2.0;
    
    vec3 color = vec3(1.0, 0.5, 0.1) * intensity;
    color += vec3(0.2, 0.4, 1.0) * pow(intensity, 2.0) * 0.2;
    
    return color;
}

vec3 aces(vec3 x) {
    const float a = 2.51;
    const float b = 0.03;
    const float c = 2.43;
    const float d = 0.59;
    const float e = 0.14;
    return clamp((x * (a * x + b)) / (x * (c * x + d) + e), 0.0, 1.0);
}

vec3 render(vec3 ro, vec3 rd) {
    vec3 pos = ro;
    vec3 vel = rd;
    vec3 accColor = vec3(0.0);
    float glow = 0.0;
    
    for (int i = 0; i < MAX_STEPS; i++) {
        float r = length(pos);
        
        if (r < 2.0 * 1.01) break;
        if (r > MAX_DIST) break;
        
        vec3 toCenter = -normalize(pos);
        float strength = 3.5 / (r * r);
        vel += toCenter * strength * STEP_SIZE;
        vel = normalize(vel);
        
        vec3 prevPos = pos;
        pos += vel * STEP_SIZE * 5.0;
        
        if (prevPos.y * pos.y < 0.0) {
            float t = abs(prevPos.y) / (abs(prevPos.y) + abs(pos.y));
            vec3 hitPos = mix(prevPos, pos, t);
            float hitR = length(hitPos);
            
            if (hitR > 2.6 && hitR < 14.0) {
                float phi = atan(hitPos.z, hitPos.x);
                vec3 diskCol = getDiskColor(hitR, phi);
                
                vec3 diskVelDir = normalize(vec3(-hitPos.z, 0.0, hitPos.x));
                float v = 0.5 / sqrt(hitR - 1.5);
                float doppler = 1.0 + dot(vel, diskVelDir) * v;
                float beaming = pow(doppler, 3.5);
                
                accColor += diskCol * beaming * 0.4;
            }
        }
        
        float distToPlane = abs(pos.y);
        if (r > 2.6 && r < 14.0) {
             glow += (1.0 / (distToPlane * 10.0 + 1.0)) * 0.002 * (1.0/r);
        }
    }
    
    accColor += vec3(0.4, 0.6, 1.0) * glow;
    
    if (length(pos) > 2.0 * 1.01) {
        vec3 skyDir = vel;
        
        float stars = pow(hash(dot(skyDir, vec3(12.989, 78.233, 45.164))), 300.0);
        accColor += vec3(stars) * 0.8;
        
        float n = fbm(skyDir * 3.0);
        float n2 = fbm(skyDir * 6.0 + vec3(u_time * 0.05));
        
        vec3 nebulaColor = vec3(0.2, 0.0, 0.4);
        nebulaColor = mix(nebulaColor, vec3(0.0, 0.3, 0.6), n);
        nebulaColor = mix(nebulaColor, vec3(0.6, 0.2, 0.1), n2 * 0.5);
        
        float intensity = n * n2 * 2.0;
        intensity = smoothstep(0.2, 0.8, intensity);
        
        accColor += nebulaColor * intensity * 0.3;
    }
    return accColor;
}

void main() {
    vec2 uv = (gl_FragCoord.xy - 0.5 * u_resolution.xy) / u_resolution.y;
    
    vec3 camPos = vec3(0.0, 1.0, -15.0);
    float rotY = (u_mouse.x / u_resolution.x - 0.5) * 6.0;
    camPos = rotateY(camPos, rotY);
    
    vec3 forward = normalize(-camPos);
    vec3 right = normalize(cross(forward, vec3(0.0, 1.0, 0.0)));
    vec3 up = cross(right, forward);
    
    vec3 dir = normalize(forward + uv.x * right + uv.y * up);
    
    vec3 col;
    vec3 dirR = normalize(forward + (uv.x + 0.003) * right + uv.y * up);
    vec3 dirG = dir;
    vec3 dirB = normalize(forward + (uv.x - 0.003) * right + uv.y * up);
    
    col.r = render(camPos, dirR).r;
    col.g = render(camPos, dirG).g;
    col.b = render(camPos, dirB).b;
    
    col = aces(col);
    
    col = pow(col, vec3(1.0/2.2));
    
    gl_FragColor = vec4(col, 1.0);
}
