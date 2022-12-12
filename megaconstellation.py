from satellite_constellation import Constellation

num_sat=22*72
num_planes=72
phase=1
inclination=53
altitude=550
ecc=0.00036
beam_width=1.3

starlink = Constellation.WalkerConstellation(num_sat,num_planes,phase,inclination,altitude,ecc,beam_width)
starlink.as_xml()