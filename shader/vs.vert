#version 330

in vec3 aPosition;
in vec3 aNormal;
in vec3 aColor;

uniform mat4 uModelViewMatrix;
uniform mat4 uNormalMatrix;
uniform mat4 uProjectionMatrix;

out vec3 vPosition;
out vec3 vNormal;
out vec3 vColor;

void main() {
    vec4 modelViewPosition = uModelViewMatrix * vec4(aPosition, 1);
    vPosition = modelViewPosition.xyz;
    vNormal = normalize((uNormalMatrix * vec4(aNormal, 1)).xyz);
    vColor = aColor;
    gl_Position = uProjectionMatrix * modelViewPosition;
}
