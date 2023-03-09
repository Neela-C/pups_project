#version 330 core

uniform vec3 lightposition;

in vec3 fragPos;
in vec3 fragColor;
in vec3 n;

out vec4 color;

void main() {
  vec3 ambient_color = vec3(0.2f, 0.2f, 0.2f);
  vec4 base_color = vec4(fragColor, 1.0);

  



  // Interpolate between normals
  vec3 n1 = normalize(n);
  vec3 n2 = normalize(cross(n, vec3(0.0, 0.0, 1.0)));
  vec3 n1_2 = normalize(mix(n1, n2, smoothstep(-1.0, 1.0, fragPos.y)));
  vec3 norm = normalize(n1_2);
  
  // first point light
  vec3 light_pos1 = lightposition;
  vec3 light_color1 = vec3(1.f, 1.f, 1.f);
  
  vec3 light_direction1 = normalize(light_pos1 - fragPos);
  float diff1 = max(dot(norm, light_direction1), 0.0);
  vec3 diffuse1 = diff1 * light_color1;

  // second point light
  vec3 light_pos2 = vec3(100.0f, 100.0f, -100.0f);
  vec3 light_color2 = vec3(1.f, 1.f, 1.f);
  vec3 light_direction2 = normalize(light_pos2 - fragPos);
  float diff2 = max(dot(norm, light_direction2), 0.0);
  vec3 diffuse2 = diff2 * light_color2;
  //// third point light
  //vec3 light_pos3 = vec3(100.0f, -100.0f, -100.0f);
  //vec3 light_color3 = vec3(0.0f, 0.0f, 0.4f);
  //vec3 light_direction3 = normalize(light_pos3 - fragPos);
  //float diff3 = max(dot(norm, light_direction3), 0.0);
  //vec3 diffuse3 = diff3 * light_color3;
  //
  //// 4th point light
  //vec3 light_pos4 = vec3(-100.0f, 100.0f, -100.0f);
  //vec3 light_color4 = vec3(0.0f, 0.0f, 0.4f);
  //vec3 light_direction4 = normalize(light_pos4 - fragPos);
  //float diff4 = max(dot(norm, light_direction4), 0.0);
  //vec3 diffuse4 = diff4 * light_color4;

  // summing the diffuse lighting contributions from both point lights
  vec3 diffuse_color = diffuse1;// + diffuse2 ;//+ diffuse3 + diffuse4;

  // combining
  color = vec4((0.1*ambient_color + 1.5*diffuse_color) * fragColor, 1.0) * base_color;
}
