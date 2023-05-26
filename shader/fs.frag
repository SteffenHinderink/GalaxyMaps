#version 330

uniform mat4 uModelViewMatrix;
uniform mat4 uNormalMatrix;
uniform mat4 uProjectionMatrix;

in vec3 vPosition;
in vec3 vNormal;
in vec3 vColor;

out vec4 outColor;

void main() {
    vec3 ca = 0.3 * vColor;
    vec3 cd = 0.7 * vColor;
    vec3 cs = 0.1 * vec3(1, 1, 1);
    float s = 10.0;
    vec3 n = normalize(vNormal);
    vec3 v = normalize(-vPosition);
    if (dot(n, v) < 0.0) {
        n = -n;
    }
    vec3 l = v;
    vec3 r = 2.0 * n * dot(n, l) - l;
    vec3 phong = (ca + cd * max(0.0, dot(l, n)) + cs * pow(max(0.0, dot(v, r)), s)) * vec3(1, 1, 1);
    outColor = vec4(phong, 1);
}
