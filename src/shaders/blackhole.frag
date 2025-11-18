#version 120

uniform vec2 u_resolution;
uniform float u_time;

// Constants
const float PI = 3.14159265359;
const int MAX_STEPS = 100;
const float MAX_DIST = 50.0;
const float STEP_SIZE = 0.1;
const float M = 1.0; // Black hole mass

// Rotate vector
vec3 rotateY(vec3 v, float angle) {
    float c = cos(angle);
    float s = sin(angle);
    return vec3(c * v.x + s * v.z, v.y, -s * v.x + c * v.z);
}

// Get background color (procedural)
vec3 getBackground(vec3 dir) {
    float theta = acos(dir.y);
    float phi = atan(dir.z, dir.x);
    
    float u = phi / (2.0 * PI);
    float v = theta / PI;
    
    // Galaxy band
    float galaxy = exp(-10.0 * (v - 0.5) * (v - 0.5));
    
    // Stars
    float stars = fract(sin(dot(dir.xy, vec2(12.9898, 78.233))) * 43758.5453);
    stars = step(0.99, stars);
    
    vec3 color = vec3(0.05, 0.05, 0.1); // Dark blue space
    color += vec3(0.8, 0.7, 0.5) * galaxy;
    color += vec3(1.0) * stars;
    
    return color;
}

void main() {
    vec2 uv = (gl_FragCoord.xy - 0.5 * u_resolution.xy) / u_resolution.y;
    
    // Camera setup
    vec3 camPos = vec3(0.0, 0.0, -15.0);
    // Rotate camera around
    camPos = rotateY(camPos, u_time * 0.5);
    
    vec3 forward = normalize(-camPos);
    vec3 right = normalize(cross(forward, vec3(0.0, 1.0, 0.0)));
    vec3 up = cross(right, forward);
    
    vec3 dir = normalize(forward + uv.x * right + uv.y * up);
    
    // Ray marching (Simplified Geodesic Integration)
    // We use a simplified step-by-step integration for real-time performance
    // in GLSL 1.20 without complex structs.
    
    vec3 pos = camPos;
    vec3 vel = dir;
    
    bool hitHorizon = false;
    
    for (int i = 0; i < MAX_STEPS; i++) {
        float r = length(pos);
        
        if (r < 2.0 * M * 1.01) {
            hitHorizon = true;
            break;
        }
        if (r > MAX_DIST) {
            break;
        }
        
        // Gravity acceleration (Newtonian approx + strong field correction term for visual bending)
        // Real geodesic equation is complex, here we approximate the bending force
        // F = -M/r^2 * (some correction)
        // Actually, let's try to be slightly more accurate with the effective potential force
        // or just use the bending formula: d(vel)/dt
        
        // 1.5 * Rs / r^3 term is a common approximation for light bending
        vec3 acc = -1.5 * (2.0 * M) * pos / pow(r, 4.0) * length(cross(vel, cross(pos, vel))); // This is heuristic
        
        // Better heuristic:
        // Light bends by 4M/r * (impact_param)
        // Let's just use a simple 1/r^2 attraction for the photon path to simulate lensing visually
        // This is NOT physics accurate but looks like it for real-time 60fps without RK4
        // For true accuracy we need the Christoffel symbols which are heavy for GLSL 1.20 loop
        
        // Let's try a simple "gravity" pull on the photon
        // Force perpendicular to velocity, towards center
        vec3 toCenter = -normalize(pos);
        float dist = r;
        
        // Bending strength
        float strength = 3.0 * M / (dist * dist); // Tuned for visual
        
        // Update velocity direction (light speed is constant 1, so we just rotate vel)
        vel += toCenter * strength * STEP_SIZE;
        vel = normalize(vel);
        
        pos += vel * STEP_SIZE * 5.0; // Move faster
    }
    
    vec3 color = vec3(0.0);
    if (hitHorizon) {
        color = vec3(0.0);
    } else {
        color = getBackground(vel);
    }
    
    gl_FragColor = vec4(color, 1.0);
}
