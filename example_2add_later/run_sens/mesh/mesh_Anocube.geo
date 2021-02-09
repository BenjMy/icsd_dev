// Gmsh project created on Fri Aug 02 10:43:00 2019
// Landfill Porto

clb=5;
cla=1;
clf=1;

Lworld = 100;
Dworld = 80;

LAno = 10;
xAno = 0;
yAno = 0;
zAno = -10;
TAno = 10;


// outer boundaries
Point(1) = {-Lworld/2, -Lworld/2, 0, clb};
Point(2) = {-Lworld/2, Lworld/2, 0, clb};
Point(3) = {Lworld/2, -Lworld/2, 0, clb};
Point(4) = {Lworld/2, Lworld/2, 0, clb};

// Line electrodes area

// boundaries of outer area
p1[]=Translate {0, 0, -Dworld} { Duplicata { Point{1}; } };
p2[]=Translate {0, 0, -Dworld} { Duplicata { Point{2}; } };
p3[]=Translate {0, 0, -Dworld} { Duplicata { Point{3}; } };
p4[]=Translate {0, 0, -Dworld} { Duplicata { Point{4}; } };
Characteristic Length {p1[0], p2[0], p3[0], p4[0]} = clb;

// Line outer area
Line(1) = {2, 1};
Line(2) = {1, 3};
Line(3) = {3, 4};
Line(4) = {4, 2};
Line(5) = {1, 5};
Line(6) = {5, 7};
Line(7) = {7, 8};
Line(8) = {8, 6};
Line(9) = {6, 5};
Line(10) = {3, 7};
Line(11) = {4, 8};
Line(12) = {2, 6};

Curve Loop(1) = {1, 5, -9, -12};
Plane Surface(1) = {1};
Curve Loop(2) = {5, 6, -10, -2};
Plane Surface(2) = {2};
Curve Loop(3) = {7, -11, -3, 10};
Plane Surface(3) = {3};
Curve Loop(4) = {11, 8, -12, -4};
Plane Surface(4) = {4};
Curve Loop(5) = {7, 8, 9, 6};
Plane Surface(5) = {5};
Curve Loop(6) = {2, 3, 4, 1};
Plane Surface(6) = {6};

// outer boundaries ANOMALY
Point(9) = {xAno-LAno/2, yAno-LAno/2, zAno-TAno/2, cla};
Point(10) = {xAno+LAno/2, yAno-LAno/2, zAno-TAno/2, cla};
Point(11) = {xAno+LAno/2, yAno+LAno/2, zAno-TAno/2, cla};
Point(12) = {xAno-LAno/2, yAno+LAno/2, zAno-TAno/2, cla};

// translate Anomaly
p1[]=Translate {0, 0, TAno} { Duplicata { Point{9}; } };
p2[]=Translate {0, 0, TAno} { Duplicata { Point{10}; } };
p3[]=Translate {0, 0, TAno} { Duplicata { Point{11}; } };
p4[]=Translate {0, 0, TAno} { Duplicata { Point{12}; } };
Characteristic Length {p1[0], p2[0], p3[0], p4[0]} = cla;


// Line Anomaly area

Line(13) = {13, 14};
Line(14) = {14, 15};
Line(15) = {15, 16};
Line(16) = {16, 13};
Line(17) = {13, 9};
Line(18) = {9, 10};
Line(19) = {10, 11};
Line(20) = {11, 12};
Line(21) = {12, 9};
Line(22) = {14, 10};
Line(23) = {15, 11};
Line(24) = {16, 12};

Curve Loop(7) = {16, 13, 14, 15};
Plane Surface(7) = {7};
Curve Loop(8) = {17, -21, -24, 16};
Plane Surface(8) = {8};
Curve Loop(9) = {13, 22, -18, -17};
Plane Surface(9) = {9};
Curve Loop(10) = {21, 18, 19, 20};
Plane Surface(10) = {10};
Curve Loop(11) = {24, -20, -23, 15};
Plane Surface(11) = {11};
Curve Loop(12) = {23, -19, -22, 14};
Plane Surface(12) = {12};

// PHYSICAL SURFACES ////// 
// 1 = Neumann No-flux, // 2= Mixed boundary conditions
// Physical Surface(2)= {1, 2, 3, 4}; // 5 = bottom of landfill area
//Physical Surface(1)= {8,9}; // 5 = bottom of landfill area


Surface Loop(1) = {5, 3, 4, 1, 6, 2};
Volume(1) = {1};
Surface Loop(2) = {12, 11, 8, 9, 7, 10};
Volume(2) = {2};

// PHYSICAL VOLUMES ////// 

//Physical Volume(1) = {volume numbers of volumes that are outside of the parameter domain};
Physical Volume(1) = {1};
Physical Volume(2) = {2};



// Add electrodes ////// 
Point	(	17	)	=	{	-35	,	-35	,	0	,	clf};
Point	(	18	)	=	{	-35	,	-30	,	0	,	clf};
Point	(	19	)	=	{	-35	,	-25	,	0	,	clf};
Point	(	20	)	=	{	-35	,	-20	,	0	,	clf};
Point	(	21	)	=	{	-35	,	-15	,	0	,	clf};
Point	(	22	)	=	{	-35	,	-10	,	0	,	clf};
Point	(	23	)	=	{	-35	,	-5	,	0	,	clf};
Point	(	24	)	=	{	-35	,	0	,	0	,	clf};
Point	(	25	)	=	{	-35	,	5	,	0	,	clf};
Point	(	26	)	=	{	-35	,	10	,	0	,	clf};
Point	(	27	)	=	{	-35	,	15	,	0	,	clf};
Point	(	28	)	=	{	-35	,	20	,	0	,	clf};
Point	(	29	)	=	{	-35	,	25	,	0	,	clf};
Point	(	30	)	=	{	-35	,	30	,	0	,	clf};
Point	(	31	)	=	{	-35	,	35	,	0	,	clf};
Point	(	32	)	=	{	-30	,	-35	,	0	,	clf};
Point	(	33	)	=	{	-30	,	-30	,	0	,	clf};
Point	(	34	)	=	{	-30	,	-25	,	0	,	clf};
Point	(	35	)	=	{	-30	,	-20	,	0	,	clf};
Point	(	36	)	=	{	-30	,	-15	,	0	,	clf};
Point	(	37	)	=	{	-30	,	-10	,	0	,	clf};
Point	(	38	)	=	{	-30	,	-5	,	0	,	clf};
Point	(	39	)	=	{	-30	,	0	,	0	,	clf};
Point	(	40	)	=	{	-30	,	5	,	0	,	clf};
Point	(	41	)	=	{	-30	,	10	,	0	,	clf};
Point	(	42	)	=	{	-30	,	15	,	0	,	clf};
Point	(	43	)	=	{	-30	,	20	,	0	,	clf};
Point	(	44	)	=	{	-30	,	25	,	0	,	clf};
Point	(	45	)	=	{	-30	,	30	,	0	,	clf};
Point	(	46	)	=	{	-30	,	35	,	0	,	clf};
Point	(	47	)	=	{	-25	,	-35	,	0	,	clf};
Point	(	48	)	=	{	-25	,	-30	,	0	,	clf};
Point	(	49	)	=	{	-25	,	-25	,	0	,	clf};
Point	(	50	)	=	{	-25	,	-20	,	0	,	clf};
Point	(	51	)	=	{	-25	,	-15	,	0	,	clf};
Point	(	52	)	=	{	-25	,	-10	,	0	,	clf};
Point	(	53	)	=	{	-25	,	-5	,	0	,	clf};
Point	(	54	)	=	{	-25	,	0	,	0	,	clf};
Point	(	55	)	=	{	-25	,	5	,	0	,	clf};
Point	(	56	)	=	{	-25	,	10	,	0	,	clf};
Point	(	57	)	=	{	-25	,	15	,	0	,	clf};
Point	(	58	)	=	{	-25	,	20	,	0	,	clf};
Point	(	59	)	=	{	-25	,	25	,	0	,	clf};
Point	(	60	)	=	{	-25	,	30	,	0	,	clf};
Point	(	61	)	=	{	-25	,	35	,	0	,	clf};
Point	(	62	)	=	{	-20	,	-35	,	0	,	clf};
Point	(	63	)	=	{	-20	,	-30	,	0	,	clf};
Point	(	64	)	=	{	-20	,	-25	,	0	,	clf};
Point	(	65	)	=	{	-20	,	-20	,	0	,	clf};
Point	(	66	)	=	{	-20	,	-15	,	0	,	clf};
Point	(	67	)	=	{	-20	,	-10	,	0	,	clf};
Point	(	68	)	=	{	-20	,	-5	,	0	,	clf};
Point	(	69	)	=	{	-20	,	0	,	0	,	clf};
Point	(	70	)	=	{	-20	,	5	,	0	,	clf};
Point	(	71	)	=	{	-20	,	10	,	0	,	clf};
Point	(	72	)	=	{	-20	,	15	,	0	,	clf};
Point	(	73	)	=	{	-20	,	20	,	0	,	clf};
Point	(	74	)	=	{	-20	,	25	,	0	,	clf};
Point	(	75	)	=	{	-20	,	30	,	0	,	clf};
Point	(	76	)	=	{	-20	,	35	,	0	,	clf};
Point	(	77	)	=	{	-15	,	-35	,	0	,	clf};
Point	(	78	)	=	{	-15	,	-30	,	0	,	clf};
Point	(	79	)	=	{	-15	,	-25	,	0	,	clf};
Point	(	80	)	=	{	-15	,	-20	,	0	,	clf};
Point	(	81	)	=	{	-15	,	-15	,	0	,	clf};
Point	(	82	)	=	{	-15	,	-10	,	0	,	clf};
Point	(	83	)	=	{	-15	,	-5	,	0	,	clf};
Point	(	84	)	=	{	-15	,	0	,	0	,	clf};
Point	(	85	)	=	{	-15	,	5	,	0	,	clf};
Point	(	86	)	=	{	-15	,	10	,	0	,	clf};
Point	(	87	)	=	{	-15	,	15	,	0	,	clf};
Point	(	88	)	=	{	-15	,	20	,	0	,	clf};
Point	(	89	)	=	{	-15	,	25	,	0	,	clf};
Point	(	90	)	=	{	-15	,	30	,	0	,	clf};
Point	(	91	)	=	{	-15	,	35	,	0	,	clf};
Point	(	92	)	=	{	-10	,	-35	,	0	,	clf};
Point	(	93	)	=	{	-10	,	-30	,	0	,	clf};
Point	(	94	)	=	{	-10	,	-25	,	0	,	clf};
Point	(	95	)	=	{	-10	,	-20	,	0	,	clf};
Point	(	96	)	=	{	-10	,	-15	,	0	,	clf};
Point	(	97	)	=	{	-10	,	-10	,	0	,	clf};
Point	(	98	)	=	{	-10	,	-5	,	0	,	clf};
Point	(	99	)	=	{	-10	,	0	,	0	,	clf};
Point	(	100	)	=	{	-10	,	5	,	0	,	clf};
Point	(	101	)	=	{	-10	,	10	,	0	,	clf};
Point	(	102	)	=	{	-10	,	15	,	0	,	clf};
Point	(	103	)	=	{	-10	,	20	,	0	,	clf};
Point	(	104	)	=	{	-10	,	25	,	0	,	clf};
Point	(	105	)	=	{	-10	,	30	,	0	,	clf};
Point	(	106	)	=	{	-10	,	35	,	0	,	clf};
Point	(	107	)	=	{	-5	,	-35	,	0	,	clf};
Point	(	108	)	=	{	-5	,	-30	,	0	,	clf};
Point	(	109	)	=	{	-5	,	-25	,	0	,	clf};
Point	(	110	)	=	{	-5	,	-20	,	0	,	clf};
Point	(	111	)	=	{	-5	,	-15	,	0	,	clf};
Point	(	112	)	=	{	-5	,	-10	,	0	,	clf};
Point	(	113	)	=	{	-5	,	-5	,	0	,	clf};
Point	(	114	)	=	{	-5	,	0	,	0	,	clf};
Point	(	115	)	=	{	-5	,	5	,	0	,	clf};
Point	(	116	)	=	{	-5	,	10	,	0	,	clf};
Point	(	117	)	=	{	-5	,	15	,	0	,	clf};
Point	(	118	)	=	{	-5	,	20	,	0	,	clf};
Point	(	119	)	=	{	-5	,	25	,	0	,	clf};
Point	(	120	)	=	{	-5	,	30	,	0	,	clf};
Point	(	121	)	=	{	-5	,	35	,	0	,	clf};
Point	(	122	)	=	{	0	,	-35	,	0	,	clf};
Point	(	123	)	=	{	0	,	-30	,	0	,	clf};
Point	(	124	)	=	{	0	,	-25	,	0	,	clf};
Point	(	125	)	=	{	0	,	-20	,	0	,	clf};
Point	(	126	)	=	{	0	,	-15	,	0	,	clf};
Point	(	127	)	=	{	0	,	-10	,	0	,	clf};
Point	(	128	)	=	{	0	,	-5	,	0	,	clf};
Point	(	129	)	=	{	0	,	0	,	0	,	clf};
Point	(	130	)	=	{	0	,	5	,	0	,	clf};
Point	(	131	)	=	{	0	,	10	,	0	,	clf};
Point	(	132	)	=	{	0	,	15	,	0	,	clf};
Point	(	133	)	=	{	0	,	20	,	0	,	clf};
Point	(	134	)	=	{	0	,	25	,	0	,	clf};
Point	(	135	)	=	{	0	,	30	,	0	,	clf};
Point	(	136	)	=	{	0	,	35	,	0	,	clf};
Point	(	137	)	=	{	5	,	-35	,	0	,	clf};
Point	(	138	)	=	{	5	,	-30	,	0	,	clf};
Point	(	139	)	=	{	5	,	-25	,	0	,	clf};
Point	(	140	)	=	{	5	,	-20	,	0	,	clf};
Point	(	141	)	=	{	5	,	-15	,	0	,	clf};
Point	(	142	)	=	{	5	,	-10	,	0	,	clf};
Point	(	143	)	=	{	5	,	-5	,	0	,	clf};
Point	(	144	)	=	{	5	,	0	,	0	,	clf};
Point	(	145	)	=	{	5	,	5	,	0	,	clf};
Point	(	146	)	=	{	5	,	10	,	0	,	clf};
Point	(	147	)	=	{	5	,	15	,	0	,	clf};
Point	(	148	)	=	{	5	,	20	,	0	,	clf};
Point	(	149	)	=	{	5	,	25	,	0	,	clf};
Point	(	150	)	=	{	5	,	30	,	0	,	clf};
Point	(	151	)	=	{	5	,	35	,	0	,	clf};
Point	(	152	)	=	{	10	,	-35	,	0	,	clf};
Point	(	153	)	=	{	10	,	-30	,	0	,	clf};
Point	(	154	)	=	{	10	,	-25	,	0	,	clf};
Point	(	155	)	=	{	10	,	-20	,	0	,	clf};
Point	(	156	)	=	{	10	,	-15	,	0	,	clf};
Point	(	157	)	=	{	10	,	-10	,	0	,	clf};
Point	(	158	)	=	{	10	,	-5	,	0	,	clf};
Point	(	159	)	=	{	10	,	0	,	0	,	clf};
Point	(	160	)	=	{	10	,	5	,	0	,	clf};
Point	(	161	)	=	{	10	,	10	,	0	,	clf};
Point	(	162	)	=	{	10	,	15	,	0	,	clf};
Point	(	163	)	=	{	10	,	20	,	0	,	clf};
Point	(	164	)	=	{	10	,	25	,	0	,	clf};
Point	(	165	)	=	{	10	,	30	,	0	,	clf};
Point	(	166	)	=	{	10	,	35	,	0	,	clf};
Point	(	167	)	=	{	15	,	-35	,	0	,	clf};
Point	(	168	)	=	{	15	,	-30	,	0	,	clf};
Point	(	169	)	=	{	15	,	-25	,	0	,	clf};
Point	(	170	)	=	{	15	,	-20	,	0	,	clf};
Point	(	171	)	=	{	15	,	-15	,	0	,	clf};
Point	(	172	)	=	{	15	,	-10	,	0	,	clf};
Point	(	173	)	=	{	15	,	-5	,	0	,	clf};
Point	(	174	)	=	{	15	,	0	,	0	,	clf};
Point	(	175	)	=	{	15	,	5	,	0	,	clf};
Point	(	176	)	=	{	15	,	10	,	0	,	clf};
Point	(	177	)	=	{	15	,	15	,	0	,	clf};
Point	(	178	)	=	{	15	,	20	,	0	,	clf};
Point	(	179	)	=	{	15	,	25	,	0	,	clf};
Point	(	180	)	=	{	15	,	30	,	0	,	clf};
Point	(	181	)	=	{	15	,	35	,	0	,	clf};
Point	(	182	)	=	{	20	,	-35	,	0	,	clf};
Point	(	183	)	=	{	20	,	-30	,	0	,	clf};
Point	(	184	)	=	{	20	,	-25	,	0	,	clf};
Point	(	185	)	=	{	20	,	-20	,	0	,	clf};
Point	(	186	)	=	{	20	,	-15	,	0	,	clf};
Point	(	187	)	=	{	20	,	-10	,	0	,	clf};
Point	(	188	)	=	{	20	,	-5	,	0	,	clf};
Point	(	189	)	=	{	20	,	0	,	0	,	clf};
Point	(	190	)	=	{	20	,	5	,	0	,	clf};
Point	(	191	)	=	{	20	,	10	,	0	,	clf};
Point	(	192	)	=	{	20	,	15	,	0	,	clf};
Point	(	193	)	=	{	20	,	20	,	0	,	clf};
Point	(	194	)	=	{	20	,	25	,	0	,	clf};
Point	(	195	)	=	{	20	,	30	,	0	,	clf};
Point	(	196	)	=	{	20	,	35	,	0	,	clf};
Point	(	197	)	=	{	25	,	-35	,	0	,	clf};
Point	(	198	)	=	{	25	,	-30	,	0	,	clf};
Point	(	199	)	=	{	25	,	-25	,	0	,	clf};
Point	(	200	)	=	{	25	,	-20	,	0	,	clf};
Point	(	201	)	=	{	25	,	-15	,	0	,	clf};
Point	(	202	)	=	{	25	,	-10	,	0	,	clf};
Point	(	203	)	=	{	25	,	-5	,	0	,	clf};
Point	(	204	)	=	{	25	,	0	,	0	,	clf};
Point	(	205	)	=	{	25	,	5	,	0	,	clf};
Point	(	206	)	=	{	25	,	10	,	0	,	clf};
Point	(	207	)	=	{	25	,	15	,	0	,	clf};
Point	(	208	)	=	{	25	,	20	,	0	,	clf};
Point	(	209	)	=	{	25	,	25	,	0	,	clf};
Point	(	210	)	=	{	25	,	30	,	0	,	clf};
Point	(	211	)	=	{	25	,	35	,	0	,	clf};
Point	(	212	)	=	{	30	,	-35	,	0	,	clf};
Point	(	213	)	=	{	30	,	-30	,	0	,	clf};
Point	(	214	)	=	{	30	,	-25	,	0	,	clf};
Point	(	215	)	=	{	30	,	-20	,	0	,	clf};
Point	(	216	)	=	{	30	,	-15	,	0	,	clf};
Point	(	217	)	=	{	30	,	-10	,	0	,	clf};
Point	(	218	)	=	{	30	,	-5	,	0	,	clf};
Point	(	219	)	=	{	30	,	0	,	0	,	clf};
Point	(	220	)	=	{	30	,	5	,	0	,	clf};
Point	(	221	)	=	{	30	,	10	,	0	,	clf};
Point	(	222	)	=	{	30	,	15	,	0	,	clf};
Point	(	223	)	=	{	30	,	20	,	0	,	clf};
Point	(	224	)	=	{	30	,	25	,	0	,	clf};
Point	(	225	)	=	{	30	,	30	,	0	,	clf};
Point	(	226	)	=	{	30	,	35	,	0	,	clf};
Point	(	227	)	=	{	35	,	-35	,	0	,	clf};
Point	(	228	)	=	{	35	,	-30	,	0	,	clf};
Point	(	229	)	=	{	35	,	-25	,	0	,	clf};
Point	(	230	)	=	{	35	,	-20	,	0	,	clf};
Point	(	231	)	=	{	35	,	-15	,	0	,	clf};
Point	(	232	)	=	{	35	,	-10	,	0	,	clf};
Point	(	233	)	=	{	35	,	-5	,	0	,	clf};
Point	(	234	)	=	{	35	,	0	,	0	,	clf};
Point	(	235	)	=	{	35	,	5	,	0	,	clf};
Point	(	236	)	=	{	35	,	10	,	0	,	clf};
Point	(	237	)	=	{	35	,	15	,	0	,	clf};
Point	(	238	)	=	{	35	,	20	,	0	,	clf};
Point	(	239	)	=	{	35	,	25	,	0	,	clf};
Point	(	240	)	=	{	35	,	30	,	0	,	clf};
Point	(	241	)	=	{	35	,	35	,	0	,	clf};
Point	(	242	)	=	{	45	,	45	,	0	,	clf};
Point	(	243	)	=	{	-45	,	-45	,	0	,	clf};
Point	(	244	)	=	{	0	,	0	,	zAno	,	clf};


// Add electrodes in surfaces
Point{17:243}	In Surface{6};	

// Add electrodes in volume
Point{244}	In Volume{2};	
//Point{244} In Surface{10};	

// PHYSICAL POINTS
// Physical Point (99) = {17:244};	





