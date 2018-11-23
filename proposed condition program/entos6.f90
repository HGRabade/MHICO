!----------------------------------------------------------------------------------------------------------------------------------------------------
!MHICO - Modelo HIdrodinámico COnjunto.
!Programa MEF válido para mallas de hasta 99999 elementos o 99999 nodos (25/02/2012). Autor: Héctor García Rábade
!Lectura de ficheros .txt (conectividades, coordenadas, propiedades y CC) y generación de ficheros .dat con formato por Tecplot	(Tecplot Data Loader)
!----------------------------------------------------------------------------------------------------------------------------------------------------

!Posibilidades del programa:
!---------------------------
!1. Se pueden resolver tres ecuaciones 2D diferentes: ec. NS2D (en función del calado en vez de la presión, ec. teórica sin conservación de masa en 3D), 
! ec. aguas someras (integración de ec. NS3D en altura, prop de manning) y ec. agua subterránea (integración de ec. continuidad en altura, propiedades 
! de conductividad, porosidad y ángulo de anisotropía) para acuíferos no confinados.
! En ec. NS2D se calculan las velocidades y el calado (igual a la altura de la lámina), en ec. aguas someras se calculan velocidades y la altura de la 
! lámina (el calado postproceso), en ec. agua subterránea se calcula la altura de la lámina o nivel freático (el espesor freático y las velocidades 
! postproceso).
!2. La ec. aguas someras incorpora un término de fricción dependiente del n=número de manning (sol. en canales con Re bajo y n=0 equivalente a sol. 
! con Re alto y cierto n).
!3. Uso de elementos de tres nodos en la ec subterránea y de elementos mixtos P2-P1 en la ec. NS2D y la ec. aguas someras.
!4. Tiene tres modelos de elementos finitos: superficial (ec. NS o aguas someras), subterráneo (ec. subt) y conjunto del tipo superficial-subterráneo (ec. 
! aguas someras y subt) conectado para acuíferos libres. 
!5. Elementos interfaz - El modelo de aguas someras utiliza la condición seco-mojado, y el modelo conjunto utiliza una condición similar que utiliza CC 
! derivadas del intercambio de información entre modelos.
!6. Modelo conjunto con conservación de masa en el contorno móvil (uso de 1 sólo contorno móvil ya que lo permite, aún en caso de trabajar con espesores 
! diferentes al haber el mismo número de nodos). 
! Cálculo para cada modelo de alturas, espesores (calado o espesor freático) y velocidades. En el modelo conjunto se representan conjuntamente 
! velocidades, alturas e individualmente el calado (valor nulo si hay espesor freático) y el espesor freático (valor nulo si hay calado). En el modelo 
! superficial se representan además la vorticidad y las tensiones.
!7. Si se simula conjunto siempre se puede simular superficial (condición seco-mojado).
!8. Posibilidad de resolver un esquema implícito (para todas las ec; resolución estacionaria conlleva a un problema no lineal, resolución transitoria 
! conlleva un problema no lineal por incremento de tiempo) o uno semi-implícito de crank-nickolson (sólo para modelo superficial, ec. NS2D o aguas 
! someras; resolución siempre transitoria, que conlleva a un problema lineal por incremento de tiempo -solución estacionaria si CC ctes-) 
!9. Posibilidad de usar el método de newton (sólo para ec. aguas someras). 
!10. Cálculo de las integrales de contorno de las tensiones viscosas (muchas veces despreciadas) para ec. NS2D o ec aguas someras
!11. Estabilización para ec. NS2D o aguas someras. Mejora la resolución del sistema lineal y la convergencia del problema no lineal.
!12. Uso de Gradientes Biconjugados con un precondicionador diagonal (usando la norma de los coeficientes de la fila, PCGB) o con un precondicionador LU 
! aproximado (PCGBLU). Cambiando el código (partes comentadas) también podría utilizarse PCGB sin precondicionamiento o con precondicionamiento con los 
! coeficientes de la diagonal.
!13. Se puede aplicar lluvia (ec. aguas someras y agua subterránea) y bombeos (ec. subt) por pantalla. 
!14. Los archivos .dat tienen el formato adecuado (en doble precisión) para que Tecplot pueda hacer una película con la solución transitoria.
!15. Para el mismo nº de elementos el modelo superficial es más lento por haber mayor nº de ecuaciones y por uso de elementos P2P1. Incremento de la 
! velocidad de computación con: 
! Utilización de vectores para el almacenamiento de la matriz dispersa del sistema (necesario en mallas de muchos elementos para tener espacio para 
! almacenar la matriz).
! Optimización en el almacenamiento de los coeficientes, mediante escritura de coeficientes en misma posición y su posterior reordenamiento de los vectores
! Uso de la última solución de Picard como aproximación para resolver el siguiente lineal.
! Modelo superficial: posible uso para Re bajos uso del método de newton (aguas someras) o uso de PCGBLU, posible uso para Re altos de estabilización.  
! Modelo subterráneo: posible uso de PCGBLU
! Modelo conjunto: uso de 1 it sup por cada sol del modelo subterráneo lo que es necesario porque permite acoplamiento y lo señalado para cada modelo 
! (si se usa un método de PCGB será para ambos modelos). 
!16. Código compatible con .f95. Ello permite su uso en compiladores de 64 bits con los que la ejecución del programa será mucho más veloz.

!Motivos de error (datos de etiquetas para fácil localización mediante búsqueda de etiquetas):
!---------------------------------------------------------------------------------------------
!1. Dependiente del ejemplo (denotado por etiqueta !ul): Dar valores de velocidad y longitud características en la fórmula del número de Re para calcular 
! Re (formulas escritas en código para cada ejemplo que fue utilizado).
!2. Dependiente del ejemplo (denotado por etiqueta !ac): Dar buena condición inicial o aproximación inicial (CI-aprox) para flujo superficial y 
! evitar definir vibv(u)=z(u) para flujo subterráneo. Dar sólo valor trivial si navier='si' o navier='no',ap='si' (altura=calado=0). De todos modos
! la selección previa evita que existan calados negativos o nulos los cuales no permitirían solución o generarían números de valor infinito (Nan). 
!3. Dependiente del ejemplo: Selección de dominio superficial a través de propiedad, no consideración de la lluvia en modelo superficial, modificación 
! de hmin (abrir o cerrar las etiquetas !es). 
!4. Mala edición de la malla: mal dadas las CC en los ficheros mallainicial.txt o mallasubinicial.txt (por ejemplo hay CC de vel=0 entre nodos con CC de 
! H=altura de la lámina en mallainicial.txt, falta de CC en algún nodo del contorno,...). En el modelo conjunto es posible prescindir de CC de nivel 
! freático para el flujo subterráneo.
!5. Dependiente del ejemplo: faltan ficheros. El fichero mallainicial es necesario para los modelos superficial y conjunto. El fichero mallasubinicial es 
! necesario para los modelos subterráneo y conjunto. El fichero propiedad.txt es por defecto necesario (a veces no se dan todas las propiedades con él). 
! En caso de no utilizarlo comentar líneas (con la etiqueta !fi). Siempre considerar dar propiedades constantes (en las líneas con la etiqueta !fa),
! propiedades no nulas con sentido físico son necesarias para el modelo subterráneo.
!6. No convergencia: La variable tol1=tol dada en la subrutina gradientesbiconjugados puede representar problemas. En la resolución del sistema lineal no 
! se debe alcanzar un error superior al que se pide en los modelos tol2=tol (del sistema no lineal). Si se alcanza se impidiría la convergencia y este error 
! suele ser mayor a tol1 (error estimado<tol1 pero error>tol1) en alguna componente. Por ello tol1 debe ser más pequeño.     
! Se ha visto que en los peores casos llega con tol1=1e-10 para tol2=1e-6 (caso programado) o tol1=1e-9 para tol2=1e-5. Además en la representación se 
! diferenciarán mejor las variaciones aunque no esté resuelto el problema no lineal con gran precisión.
!7. Evitable en ejecución (usuario elige por pantalla): Se están usando las integrales de contorno de términos viscosos (también al usar el método 
! de Newton) y esto da problemas.
!8. Evitable en ejecución (usuario elige por pantalla): Si no se logra solución para el sistema lineal, el precondicionador para resolverlo puede no 
! ser adecuado y se aconseja elegir el otro (a través de PCGB - subrutina gradientesbiconjugados ó PCGBLU - subrutina dslubc, será más rápido el LU). 
! En el modelo subterráneo cualquier precondicionador funcionará bien. En el superficial: el buen funcionamiento de cada precondicionador dependerá del 
! ejemplo. 
! No se da la opción del uso del precondicionador LU en caso de usar el método de Newton o la estabilización (modificar el código para ello) porque se 
! han visto peores resultados.
!9. Improbable. El precondicionador diagonal podría no funcionar bien en el modelo superficial porque podrían aparecer números infinitos (Nan en la 
! variable di). De todos modos, está programado para que no suceda (ver sub nuevamalla en el caso particular de considerar ciertos elementos).  
!10. No convergencia: En caso de que el modelo superficial no converja (para Re altos pasa) considerar usar estabilización (sistema mejor condicionado). 
! Ésta es adecuada si los elementos de la malla tienen forma regular (misma forma, puediendo tener diferente tamaño) y funcionará mejor tanto más similar 
! sea la forma.
!11. No convergencia del método de Newton: Si se está usando el método de newton y no converge tal vez llegue con subir el número de iteraciones previas 
! de Picard en el código. 
!12. Imposible. Al usar el esquema estabilizado en transitorio, los términos de masa con pesos SUPG-PSPG generan problemas, sobre todo si se dan incrementos 
! de tiempo pequeños (los términos van dividos por estos incrementos). Cambiando el código (partes comentadas) se pueden utilizar.
!13. No convergencia inducida: hmin=1e-3 podría no ser suficientemente grande como para evitar problemas de localización y deslocalización de los mismos 
! elementos indefinidamente. 
!14. Soluciones con error: Si se aplica el modelo superficial con un esquema semi-implícito utilizar un incremento de tiempo pequeño o habrá errores en la 
! solución. Tampoco se recomienda usar un valor muy pequeño porque tardaría demasiado (aún no se ha limitado el incremento de tiempo con el parámetro CFL).

!Notas sobre la aprox-CI y la selección de dominios:
!---------------------------------------------------
!1. El fichero mallainicial.txt lleva para el contorno de todo el dominio unas CC donde hay flujo superficial y CC nulas donde no lo hay (CC superficiales). 
! El fichero mallasubinicial.txt lleva para el contorno de todo el dominio unas CC donde hay flujo subterráneo y CC nulas donde no lo hay (CC subterráneas). 
!2. NS2D 
! La ecuación no tiene en cuenta cota de fondo variable y el proceso iterativo permite calados negativos.
! Está dada una aproximación inicial con solución trivial para estacionario, le vale cualquiera.
! El transitorio necesita de una CI=vib (uso de vib) y también vale la solución trivial. 
!3. Aguas someras 
! Poner aproximación para arrancar el estacionario (uso de vib). 
! Con ella, además se efectúa una selección de dominio pues el proceso iterativo no permite calados negativos (luego se hace cargo la condición seco-mojado).
! Otra opción es usar una solución tras una (o más) iteración de la ec NS2D (ap='si'). Aquí no es necesaria una aprox específica (se usa vib) y se efectuará
! la selección del dominio del mismo modo.
! En caso donde la lámina de agua pueda ser prácticamente horizontal (CC(H)>zmax de zona seleccionada) vale cualquier opción. 
! Es razonable una aprox con altura contante y velocidades nulas que provoque calados. 
! En otro caso, para la primera opción se hace necesario tener una solución real (sol con calados y velocidades, lo que es imposible) por lo que se usará 
! la segunda opción con un Re que permita calados en un dominio seleccionado. Se recomienda seleccionar el dominio a través de una propiedad tipo 
! coef de manning (es casi imposible que una aprox-CI inventada permita una selección realista).
! Truco: subir hmin en nuevamalla (c. seco-mojado), hacer más it previas de NS (por ejemplo hasta que no haya variación de elementos de dominio superficial 
! entre 2 it), buscar con N-S 2D un calado constante (con NS2D las velocidades son conservativas en 2D indep del calado y para un calado cte lo serán en 3D),
! cuidar que no se seque la zona principal con flujo por c. seco-mojado, usar un coef. Manning pequeño (o nulo).  
! Poner una CI para arrancar el transitorio (uso de vib). Sólo se podrá utilizar N-S 2D (ap='si') en el primer incremento (se usa vib).
!4. Aguas someras con uso de aprox inicial con altura constante (casos donde no sea necesaria una selección a través de una propiedad): si se resuelve 
! el caso estacionario con el modelo superficial o el conjunto (con aguas someras), y no se quiere resolver inicialmente en todo el dominio para las 
! CC superficiales se necesita una aprox inicial con menor cota para la altura de agua. 

!Modelo conjunto (para cada incremento de tiempo o para el transitorio se empieza por el modelo superficial):
!------------------------------------------------------------------------------------------------------------
!1. Condiciones de contorno: Se deben diferenciar las CC de flujo superficial de las CC de vel subterráneas aplicadas como CC superficial.
! Para simular flujo superficial se ve la necesidad de imponer velocidades superficiales y calados. 
! Las citadas velocidades subterráneas apenas generarán flujo superficial (no se obtiene el flujo propio de un río con CC de calado y 
! condiciones de este tipo, de infiltración). Con ellas es posible simular zonas de agua embalsada (flujo superficial aislado), habiendo 
! CC de vel subterránea únicamente en el contorno móvil. Además la solución tendrá sentido físico (si sólo hay infiltración).
!2. Flujo aislado: El modelo no genera zonas aisladas superficiales nuevas, sólo tiene en cuenta las zonas aisladas superficiales localizadas 
! previamente o las localizadas dentro del subdominio superficial (por ejemplo si en una simulación transitoria la altura de agua decrece). 
! El flujo superficial puede generar subdominios subterráneos (flujo subterráneo aislado) durante el proceso iterativo.
! El flujo subterráneo no puede generar subdominios superficiales (flujo superficial aislado) durante el proceso iterativo.
! Las zonas de flujo superficial aislado se localizarán con la CI-aprox o una propiedad (p.e. número de manning) y si se deslocalizan no se vuelven 
! a tener en cuenta (la consideración de todo el dominio por la CI-aprox obligará a que las alturas se ajusten a la altura definida como CC superficial).
! Si existe una zona así, ahí el agua debe encontrarse contenida por una depresión del terreno para asegurar que no haya vertido 
! (seleccionar una buena CI-aprox). Si esto no se cumple (puede ocurrir para un tiempo de simulación determinado), la zona aislada se acabará uniendo 
! a otra zona de forma errónea ya que la condición seco-mojado sumará elementos sin considerar el incremento de tiempo. Además no se tendrán condiciones 
! de flujo superficial en el contorno de la zona, el agua no tendrá la suficiente velocidad como para que i sea más o menos I existiendo 
! calado normal, por lo que no se obtendrá una solución del tipo de un afluente bajando a un río.
!3. Se ha visto que el modelo conjunto no permite el cálculo de una solución estacionaria si se está simulando flujo superficial aislado (sin CC de flujo 
! superficial en el contorno de todo el dominio). Tal vez porque la solución estacionaria conlleve al vertido (y sean necesarias esas CC).
!4. La exactitud en la velocidad subterránea impuesta en el contorno móvil dependerá del tamaño de los elementos interfaz pegados al contorno ya que
! es calculada a través de un nivel freático medio (un plano). El caudal que se intercambia entre los modelos es función de ella.
 
!Posibilidades (etiquetas !cñ):
!------------------------------
!1. Existen dos posibilidades para dar CC (intersección contorno móvil-contorno) en el modelo conjunto desde la subrutina mallaaguassomeras.
!2. Postproceso existen dos formas de calcular el calado en nodos cuadráticos pertenecientes al dominio subterráneo cercanos al contorno móvil.
!3. Posibilidad de calcular la solución subterránea con una matriz de masa concentrada (teóricamente mejora el sistema y no tiene error en el 
!caso estacionario).
!4. Con la variable yn se podría calcular una solución para aguas someras con ix (pend en x)=Ix e iy (pend en y)=Iy independiente de la velocidad. 
! Siempre generará soluciones con el mismo calado en todos los nodos independientemente de la cota del terreno con cierto sentido en el caso de 
! canales con igual cota de fondo en cada sección transversal (canal rectangular). En este caso se generaría un calado normal para cualquier velocidad. 
!5. Cálculo preproceso de pendientes medias (a través de los diferentes valores nodales que se tienen) en direcciones x e y multiplicadas por una 
! constante. No utilizada.

!Desarrollos futuros:
!--------------------
!1. Programar cogiendo CI de fichero y CC variables (de momento se utilizan CC constantes en el tiempo). En caso de usar CC variables 
! habría que sustituir las CC correspondientes a ese instante en la solución en el instante anterior (por ejemplo quitando la condición de itt=0 
! en la subrutina aguassomeras)	para el cálculo de las integrales de contorno.
!2. Considerar la lluvia en el modelo de aguas someras estabilizado.
!3. Cálculo de más variables post-proceso. Algunas posibilidades están comentadas (etiqueta !ci).
!4. Considerar aplicar el método de Newton a las ecuaciones N-S 2D (como en Reddy y Gartling).
!5. Resolver la ecuación de convección-difusión para el cálculo del transporte de solutos (utilizando los resultados hidrodinámicos). 
!6. En subrutina aguassubterranea sacar del bucle (que aplica el método de Picard) el cálculo de vnin, dimvectsb y las cajas lineales del tiempo 
! (no necesario al no variar el dominio).
!7. No hay posibilidad de definir la intensidad de lluvia por regiones. Incorporar la lluvia y bombeo en el fichero de propiedades.
!8. Comprobar si no es necesario ordenar coeficientes en la subrutina dsilus.  

!-------------------------------------------------------------------------------------------------------------------------------------------------------
!INICIO DE PROGRAMA
!-------------------------------------------------------------------------------------------------------------------------------------------------------

!Módulo para los vectores de mayor dimensión
!---------------------------------------------
!Los módulos permitirán generar variables globales que pueden ser utilizadas y modificadas desde cualquier subrutina.
module allocatacion
!Al igual que en ciertas subrutinas, se utiliza dimensionamiento dinámico (se dará la dimensión en tiempo de ejecución con allocate).
!Las siguientes serán las variables que usan más memoria. Se han usado las mínimas posibles.
integer*4, dimension(:),allocatable::ita,cia,cja 			 
real*8, dimension(:),allocatable::sa,ca
endmodule

!Módulo para las subrutinas donde se produce intercambio de condiciones de contorno
!------------------------------------------------------------------------------------
module interaccion
integer*4 b,c,e,f,sb(4),eaf,no(6)								  
character a*80
endmodule

!Módulo para las subrutinas donde se calculan las matrices elementales
!-----------------------------------------------------------------------
module elemental
integer*4, dimension(:),allocatable::pos,posi,posdosi,postresi
integer*4 u,uu,uuu,ui,uj,n(6)								  
real*8 xa,ya,xb,yb,luno(7),ldos(7),ja(7),lunot(7),ldost(7),jat(7),jac
character ac*80
endmodule

!------------------------------------------------------------------------------------------
!Programa principal
!'i' es el nº de nodos, 'j' el nº de elementos, 'modelo' gestiona que modelo se resuelve	 
!------------------------------------------------------------------------------------------
program flujo2D
!Se usa un módulo como se hará en otras subrutinas.
use elemental
!Se definen variables globales (programa principal).
integer*4, dimension(:),allocatable::vn 
integer*4 i,j 
character modelo*12	

12  format(/,I5)
23  format(4/,A80)
26  format(6X,6(X,I5))

!Fichero donde se escribirá la solución	y fichero donde se escribirán las propiedades.
open(unit=8,file='C:\solucionaguassomeras-subterraneo.dat',status='unknown')
open(unit=10,file='C:\solucionpropiedades.dat',status='unknown')	 

!Selección del modelo a utilizar.
!En caso de ser necesaria la condición seco-mojado para el modelo superficial o la condición similar para el modelo conjunto
!se trabajará con mallas menores a la malla de todo el dominio (en los ficheros malla.txt y mallasub.txt) para el cálculo de los coeficientes
!generando sólo las matrices elementales necesarias. En este caso el programa trabajará siempre con la misma numeración de nodos y se darán CC nulas 
!fuera de la malla para construir un sistema menor.
write(6,*)'Escribe el modelo a usar: superficial/subterraneo/conjunto'
read(5,*) modelo
if (modelo.eq.'superficial') then
 write(8,'(A)')'TITLE = "Aguassomeras"'
 write(8,'(A)')'VARIABLES =  X, Y, Htzt, H, VELX, VELY, VEL, TENX, TENY,'
 write(8,'(A)')'TENXY, VORT' 
 write(6,*)'Debe estar en C:\ el fichero mallainicial.txt'	 
 open(unit=1,file='C:\mallainicial.txt',status='old')
elseif (modelo.eq.'subterraneo') then
 write(8,'(A)')'TITLE = "Subterraneo"'
 write(8,'(A)')'VARIABLES =  X, Y, Hd, E, VELX, VELY, VEL'	  
 write(6,*)'Debe estar en C:\ los ficheros mallasubinicial.txt'	 
 open(unit=1,file='C:\mallasubinicial.txt',status='old')
else
 write(8,'(A)')'TITLE = "Aguassomeras y subterráneo"'
 write(8,'(A)')'VARIABLES =  X, Y, HtHd, Htzt, H, E, VELX, VELY, VEL' 
 write(6,*)'Deben estar en C:\ los ficheros mallainicial.txt, mallasubinicial.txt'	
 open(unit=1,file='C:\mallainicial.txt',status='old')
endif

!Lectura de número de nodos y elementos de la malla desde el fichero existente.
read (1,12) i
read (1,12) j
rewind(1)

!En los vectores pos se guardará la posición donde empieza la siguiente caja para cada fila.
allocate (vn(i),pos(i),posi(3*i),posdosi(3*i),postresi(3*i))	 

do u=1,i
vn(u)=0
enddo
read(1,23)ac 
 do u=1,j
  !Si vn es 1 el nodo es cuadrático. Así vn llevan los nodos "no esquina" de toda la malla.
  read(1,26)n(1),n(4),n(2),n(5),n(3),n(6)	  
  do uu=1,3
  vn(n(3+uu))=1
  enddo
 enddo 
close(1)
	  
call aguassomerassubt (i,j,vn,modelo)
deallocate(vn,pos,posi,posdosi,postresi)
close(8)
close(10) 
end

!-----------------------------------------------------------------------------------------------------------------------------------------------------------
!Subtrutina AGUASSOMERASSUBT donde se calcula la solución para flujo superficial, subterráneo o conjunto.
!'i' es el nº de nodos, 'j' el nº de elementos, l el número de nodos esquina (de la toda la malla), 
!'ma' el número de Manning, 'nu' la viscosidad (leídas de fichero en 'zzn'), 're' el número de Reynolds, 'x,y' las coordenadas 2D de cada nodo, 
!'z' la cota del terreno en cada nodo, 'zp' la cota del sustrato impermeable en cada nodo, 'zzp' esta cota pero modificada en el contorno móvil, 
!'kx,ky,ag,nd' las permeabilidades en las direcciones x e y, el ángulo de anisotropía y la porosidad del terreno (leídas de fichero en 'zz,zzz,zza y zzzz'),
!'qx qy' son los caudales para el flujo subterráneo, 'velx,vely' las velocidades del flujo subterráneo (obtenidas postproceso),    
!'vn,ve,ven,vo,vnv,vov' son contadores de nodos cuadráticos, lineales o con condición de contorno ('vnv,vov' para flujo subterráneo), 
!'vt,vb,vtv,vbv' son los vectores donde se guarda la solución del sistema ('vtv,vbv' para flujo subterráneo),
!'vit,vib,vitv,vibv' son los vectores solución obtenidos en el instante de tiempo anterior ('vitv,vibv' para flujo subterráneo),
!'sino' permite la elección del caso estacionario o no estacionario,
!'At,Ata,tac' llevan el valor del incremento de tiempo, el valor de tiempo simulado antes de imponer más tiempo y el valor total de tiempo que se simulará,
!'tiempo,nt' el tiempo de simulación y el número de incrementos de tiempo que supone,
!'newton,bcg' permite seleccionar los métodos de Netwon/Picard, el solver BCG con precondicionador diagonal o LU,
!'navier,est,imp' permite seleccionar las ecuaciones N-S 2D o Aguas someras (modelo superficial), formulación BG o SUPG-PSPG, implícito o semi-implícito,
!'ap' permite usar la ecuación de N-S 2D en las primeras iteraciones cuando se usan las ecuaciones de aguas someras,
!'ten' permite el cálculo de las integrales de contorno de los términos viscosos, 'modelo' gestiona que modelo se resuelve, 
!'vib,vit,vibv, vitv' guardan los valores de la iteración anterior en el tiempo (vibv,vitv para flujo subterráneo),
!'qb' lleva las condiciones de contorno de bombeo para flujo subterráneo, 
!'ql' lleva las condiciones de contorno de lluvia para flujo subterráneo y superficial.
!'te,vor' la tensión y vorticidad calculada postproceso (sólo con los resultados del modelo superficial).
!'it,itt' llevan el nº de iteración del mét. para resolver la no linealidad de las ec. flujo superficial (Newton o Picard) y el nº de la iteración temporal,  
!'u,ui,ut' contador para bucles generales y contador para bucle temporal. 
!'frec,ni' permite la impresión de resultados cada cierto número de incrementos de tiempo. 
!-----------------------------------------------------------------------------------------------------------------------------------------------------------
subroutine aguassomerassubt (i,j,vn,modelo)
use interaccion
!Se definen nuevas variables, que serán locales (no traídas desde otra subrutina o el programa principal). Las variables locales deberían
!desaparecer (de la memoria) al terminar una subrutina (aunque aquí ya se hace a través del dimensionamiento dinámico para los vectores). 
integer*4 i,j,u,ui,ut,vn(i),it,itt,ss,l,nt,ni,frec				 		  	 
real*8, dimension(:),allocatable::x,y,z,zp,zzp,vt,vit,vib,vb,vo,eval,vov,vtv,vitv,vibv
real*8, dimension(:),allocatable::kix,kiy,ag,nd,ma,qx,qy,qb,ql,velx,vely,vx,vy,ht,evol,te,vor  
real*8 zu,tol,re,nu,mu,At,Ata,tiempo,tac,zz,zzz,zza,zzzz,zzn,qbb,qll  
character modelo*12,sino*2,newton*2,bcg*2,ten*2,Atee*10,Atai*29,ww*1,navier*2,est*2,imp*2,ap*2
!Se utiliza dimensionamiento dinámico aunque no sea necesario (se podría definir como la variable traída vn, que además es global). 
!Sin embargo es el modo correcto de proceder para para que la subrutina pre-asigne la memoria necesaria para las nuevas variables locales.
allocate(x(i),y(i),z(i),zp(i),zzp(i),vt(3*i),vit(3*i),vib(3*i),vb(3*i),vo(3*i),eval(i),vov(i),vtv(i),vitv(i),vibv(i))
allocate(kix(i),kiy(i),ag(i),nd(i),ma(i),qx(i),qy(i),qb(i),ql(i),velx(i),vely(i),vx(i),vy(i),ht(i),evol(i),te(3*i),vor(i))
	  			 											   
!Formatos de lectura o escritura:
!--------------------------------
15  format(6X,3(2X,F12.4))			 
17  format(X,I5,5(4X,F14.12))
18  format(X,I5,X,A1,X,F11.7)		  
21  format(3/,I5)
!Se escribirán el nº de elementos equivalentes (habrá un máximo de 399996 elementos P2-P1, por lo que llega con 6 cifras), de tres nodos.
22  format('ZONE T=',A29) 
23  format('I=',4x,I5,',',2x,'J=',4x,I6,',',x,'    DATAPACKING = POINT, ZONETYPE = FETRIANGLE')	  
24  format(2x,'ZONE  I=',4x,I5,',',2x,'J=',4x,I6,',',x,'    DATAPACKING = POINT, ZONETYPE = FETRIANGLE')
!Así, se escribirán tres columnas con a b c (siendo a,b,c los tres nodos de cada elemento).
!También vale: ...,x,'    F=FEPOINT') pero luego hay que escribir cuatro columnas con a b c c cambiando el formato 41 a format(4(3x,I5)).
!Necesario representarla con FEPOINT si se representa en Tecplot 9.0. En versiones siguientes es mejor hacerlo con FETRIANGLE (ahorra espacio de fichero). 																											  
26  format(6X,6(X,I5))
27  format(4/,A80)
!Para el fichero de resultados:
!Se escribirán las coordenadas con el mismo formato con el que se leen, con 7 enteros (o 6 y signo negativo) y 4 decimales.
!Se escribirán las propiedades con 4 enteros(o 3 y signo negativo) y 10 decimales.
38  format(2(x,F12.4),5(2x,F15.10))  
39  format(2(x,F12.4),5(2x,F15.10))   
40  format(2(x,F12.4),7(2x,F15.10))  
41  format(3(3x,I5))
48  format(3(x,F12.4),(2x,F15.10))	  
49  format(3(x,F12.4),3(2x,F15.10))  
50  format(4(x,F12.4),4(2x,F15.10))  
	  
!Inicialización de variables:
!----------------------------
!Se utiliza doble precisión para los cálculos y es necesario escribir los números reales de forma que se considere
!doble precisión. Así, se utiliza la notación 0.1125_8, que es equivalente a 1.125*1d-1, para escribir el número 0.11250000000...		
itt=0
Ata=0.0_8
tac=0.0_8
do u=1,3*i
vt(u)=0.0_8
vb(u)=0.0_8
enddo
do u=1,i
z(u)=0.0_8
zp(u)=0.0_8
zzp(u)=0.0_8 
vtv(u)=0.1_8
velx(u)=0.0_8
vely(u)=0.0_8
eval(u)=0.0_8
ql(u)=0.0_8
qb(u)=0.0_8
!Valores para la primera iteración del método para resolver la no linealidad.
!Los siguientes valores sólo se usan en caso de usar el modelo superficial con N-S 2D (valen valores nulos) o el modelo subterráneo (no valen).
!En este caso, al buscar sol. estacionaria o transitoria se usa la solucion trivial (vt=0) para N-S 2D. 
!En otro caso, al buscar sol. estacionaria o transitoria, se escribe después el valor de la aproximación inicial 
!o la CI (estará en vib), tal que vit=vib-z/vit=vib,vb=vib,vitv=vibv.
!Valores de velocidad para flujo superficial.
vit(u)=0.0_8
vit(i+u)=0.0_8
!Valores de espesor (flujo superficial) y altura (flujo subterráneo).
vit(2*i+u)=0.0_8			   
vitv(u)=0.1_8
enddo
!Se podría usar también el comando data para 'sb' con el que se definen datos que no cambian. No es el caso de b,c,e y por tanto 
!no se podría definir data b/0/,c/0/,e/0/ (equivalente a: data b,c,e/0,0,0/). 
sb=(/0,0,3,3/)
b=0
c=0
e=0
f=0
nt=0
ni=1

!Selección de estabilización, tipo de esquema y ecuaciones superficiales a utilizar.
!-----------------------------------------------------------------------------------
!Si se resuelve con el modelo subterráneo siempre se utiliza un esquema implícito.
if (modelo.eq.'superficial')then
 write(6,*)'modelo navier stokes en vez de aguas someras? (si/no)'
 read(5,*)navier
 !Resolución con/sin estabilización para aguas someras y para Navier-Stokes 2D.
 est='si'
 if (est.eq.'si') then
 write(6,*)'Ojo, resolucion con estabilizacion. No aplicar lluvia si SW (no programado)'
 endif
 !Resolución con esquema implícito/semi-implicito (incondicionalmente estables) para aguas someras y para Navier-Stokes 2D.
 !Cálculo de solución estacionaria: se llegará a las mismas soluciones si CC constantes	y es posible una solución estacionaria.
 !Con implícito (resolución transitoria) o semi-implícito se puede calcular con incrementos de tiempo (también si fuese explícito).
 !En implícito también se puede calcular estacionario ignorando términos temporales (resolución estacionaria).   
 !Cálculo de solución transitoria: En ambos casos se requiere solución con incrementos de tiempo.
 imp='si'
 if (imp.eq.'no') then
 write(6,*)'Ojo, resolucion con semi-implicito'
 endif
elseif (modelo.eq.'conjunto')then
 !Siempre se usará aguas someras.
 navier='no'
 !Resolución con/sin estabilización para aguas someras.
 est='si'
 if (est.eq.'si') then
 write(6,*)'Ojo, resolucion con estabilizacion. No aplicar lluvia en SW (no programado)'
 endif
 !Siempre se usará un esquema implícito.
 imp='si'
endif 
ap='no'

!Lectura de propiedades:
!------------------------
!Cálculo de la variable l
l=i-count(vn.eq.1)
!fi (hasta close(12) inclusive)
eaf=0
open(unit=12,file='C:\propiedad.txt',status='old',iostat=eaf)
if (eaf.eq.0)then
read(12,'(A)')a    
 !Valores de conductividad en m/s, de ángulo de anisotropía en radianes, de porosidad (adimensional) y de manning.
 !Lectura de propiedades desde fichero (conductividad y porosidad necesarias para flujo subterráneo) 
 do u=1,l	   
  read(12,17)ss,zz,zzz,zza,zzzz,zzn
  kix(ss)=zz
  kiy(ss)=zzz
  ag(ss)=zza
  nd(ss)=zzzz
  ma(ss)=zzn
 enddo
else				 
write(6,*)'Ojo, falta propiedad.txt modificar codigo segun sus especificaciones de error'
endif
close(12)

!Lectura de coordenadas y condiciones de contorno para el modelo superficial:
!----------------------------------------------------------------------------
!Línea a línea del fichero, empieza a leer las coordenadas x,y,z de cada nodo 
!en la línea nºelementos+6 del fichero para flujo superficial.
if ((modelo.eq.'superficial').or.(modelo.eq.'conjunto')) then 
open(unit=1,file='C:\mallainicial.txt',status='old')
do u=1,6+j
read(1,'(A)')a
enddo
 !Valores de coordenadas en metros (incluyendo la cota del terreno).
 do u=1,i
 read(1,15) x(u),y(u),z(u) 
 !fa (hasta ma(u) inclusive)	
 !Si no existe fichero de propiedades no se da propiedad alguna, y se debe dar un valor para Manning (vale dar valor nulo) para inicializar la variable.
 !Si existe puede no haberse dado esta propiedad (si todos los valores leídos son nulos).
 !Si existe y se tienen valores para la propiedad aquí se pueden sobreescribir. 
 !Un valor constante con sentido físico (para hormigón):
 !ma(u)=0.015_8 
 enddo
read(1,'(A)')a
do u=1,3*i
vo(u)=sqrt(2.0_8)
enddo
!Lectura de las condiciones de contorno para flujo superficial desde fichero.
!El programa lee ficheros donde hay "i" líneas con condición para la velocidad en x, e "i" líneas en y. Se permite cualquier nº líneas con condición 
!para el calado (por ejemplo tantas como nodos esquina como se tiene por defecto). 
!Se obtiene el vector de valores conocidos para cada iteración (CC Dirichlet constantes).
u=0
eaf=0
do while (eaf.eq.0)
 u=u+1	   
 read(1,18,iostat=eaf)ss,ww,zu
 if (eaf.eq.0) then			
  if (ww.eq.'o') then
    cycle
  elseif (u.le.i) then
    vo(ss)=zu
  elseif ((u.gt.i).and.(u.le.(2*i))) then
    vo(i+ss)=zu
  else
    vo(2*i+ss)=zu+z(ss)
  endif
 endif 
enddo
!Se guarda el fichero con toda la malla 'mallainicial.txt' en otro archivo (para no modificar los ficheros de datos originales si se generan
!mallas de menos elementos). El nuevo archivo será el que se utilize durante la simulación.
rewind(1)
open(unit=4,file='C:\malla.txt',status='unknown')
 do u=1,5+j
 read(1,'(A)')a 
 write(4,'(A)')a  
 enddo
close(4)
close(1)
!Valor del número de Reynolds dado por pantalla. 
write(6,*)'Introduce un valor para el numero de Reynolds: '
read(5,*) re
!Calculo de la viscosidad para ese número de Reynolds (nu=vel*L/re). Depende del ejemplo.
!'vel' es la velocidad característica (para calcular la viscosidad). Será el valor de la velocidad promedio condicion de contorno. 
!'L' es la longitud característica. Será el ancho en la dirección perpendicular a la dirección de avance del fluido. 
!En ejemplos utilizados: nu=(1.28_8*30.0_8)/re  nu=(0.256_8*0.43_8)/re  nu=(1.0_8*400.0_8)/re  nu=(0.666666*2.0_8)/re  nu=(0.33_8*6.0_8)/re   !ul
nu=(0.33_8*6.0_8)/re	
mu=nu
 !Resolución con Newton/Picard para las ecuaciones de aguas someras con esquema implícito y formulación BG. 
 if ((est.eq.'no').and.(navier.eq.'no').and.(imp.eq.'si').and.(modelo.ne.'subterraneo')) then
 write(6,*)'Quieres aplicar el metodo de newton para el modelo superficial?(si/no)'
 read(5,*)newton
 else
 newton='no'
 endif
endif

!Lectura de coordenadas y condiciones de contorno para el modelo subterráneo:
!----------------------------------------------------------------------------
!Línea a línea del fichero, empieza a leer las coordenadas x,y,zp de cada nodo 
!en la línea nºelementos+6 del fichero para flujo subuterráneo.
if ((modelo.eq.'subterraneo').or.(modelo.eq.'conjunto')) then         
open(unit=3,file='C:\mallasubinicial.txt',status='old')
do u=1,6+j
read(3,'(A)')a
enddo
 !Valores de coordenadas en metros (incluyendo cota del sustrato)
 do u=1,i
  read(3,15) x(u),y(u),zp(u) 
  !fa  (hasta nd(ss) inclusive)
  !Si no existe fichero de propiedades no se da propiedad alguna, y se deben dar valores de kix,kiy,ag y nd para inicializar las variables.
  !Si existe puede no haberse dado alguna propiedad (todos los valores leídos son nulos).
  !Si existe y se tiene una propiedad aquí se pueden sobreescribir. 
  !Propiedades de valor constante con sentido físico (p.e.):
  !kix(ss)=0.00001_8
  !kiy(ss)=0.00001_8
  !ag(ss)=0.0_8
  !nd(ss)=0.01_8
 enddo 
read(3,'(A)')a
do u=1,i
vov(u)=sqrt(2.0_8)
qx(u)=sqrt(3.0_8)					 
qy(u)=sqrt(3.0_8)
enddo
u=0
eaf=0
!Lectura de las condiciones de contorno desde fichero.
!En el fichero han de aparecer "i" líneas para las CC de caudal por metro lineal aunque sólo se tendrán valores para los nodos esquina.
do while (eaf.eq.0)
u=u+1
 read(3,18,iostat=eaf)ss,ww,zz
 if (eaf.eq.0) then			   
  if (ww.eq.'o') then
  cycle 
  elseif (u.le.i) then
    qx(ss)=zz
  elseif ((u.gt.i).and.(u.le.2*i)) then
    qy(ss)=zz
  else
    vov(ss)=zz
  endif
 endif 
enddo 						  
!Igual que antes, se guarda el fichero con toda la malla 'mallasubinicial.txt' en otro archivo.
rewind(3)
open(unit=7,file='C:\mallasub.txt',status='unknown')
 do u=1,5+j
 read(3,'(A)')a
 write(7,'(A)')a
 enddo
close(7)
close(3)
!Se imponen condiciones de bombeo desde pantalla (constantes en el tiempo).
write(6,*)'Existe algun bombeo?(si/no)'
read(5,*)sino
if (sino.eq.'si') then
write(6,*)'Introduce nodo y valor de caudal:'
read(5,*) u,qbb 
qb(u)=qbb
endif
endif

!Selección del precondicionador, y del uso de las tensiones viscosas:
!-------------------------------------------------------------------- 
if ((newton.eq.'no').and.(est.eq.'no')) then
!Sólo en los casos donde funciona bien se permite precondicionador LU. Así, sólo si se usa formulación BG y Picard se entraría, lo que 
!siempre ocurrirá para el modelo subterráneo (aunque est no tiene valor en este caso).
write(6,*) 'Utilizar PBCGLU en vez de PCGB (si/no)?'
read(5,*)bcg
else
bcg='no'
endif
if ((modelo.eq.'conjunto').or.(modelo.eq.'superficial')) then
write(6,*) 'Utilizar tensiones viscosas (modelo superficial) (si/no)?'
read(5,*)ten
endif 

!Condiciones iniciales o aproximación inicial   !ac
!--------------------------------------------------
!En caso de usar el modelo superficial con N-S 2D o el modelo subterráneo, sólo será necesario definir vib ó vibv como CI si se simula en transitorio 
!(se usarán valores diferentes de CI y de valor para arrancar el método de resolución del sistema no-lineal, a no ser que se escriban los mismos). 
!En este caso habrá aprox inicial válida en vit para el estacionario (para arrancar el método). Es así porque no se hará vit=vib, vitv=vibv. 
!En otro caso, siempre será necesario escribir aquí el valor de vib ó vibv. En transitorio será la CI (mismo valor de CI y de valor para arrancar el 
!método de resolución del sistema no-lineal). En estacionario será la aproximación inicial (para arrancar el método). Es así porque se 
!hará vit=vib, vitv=vibv. 
!Flujo superficial: la iteración no lineal requerirá solución con calados (vit) y la iteración temporal requerirá solución con alturas (vib).
do u=1,i		        
!Valores de velocidad para flujo superficial.
vib(u)=0.0_8           
vib(i+u)=0.0_8         
!Valores de altura de la lámina de agua para flujo superficial (en vib) y de altura de agua o nivel freático para flujo subterráneo (vibv).
!Se debería dar el mismo valor inicial para ambos modelos de forma que vib(2*i+u)=vibv(u). 
!Ojo con definir también vibv si se usa el modelo superficial de aguas someras por la posterior selección de dominios (es recomendable 
!definir vibv=vib si se desea utilizar todo el dominio).
!Ojo si vibv(u)=z(u) porque con la posterior selección de dominios aquellas zonas que no pertenezcan cumplirán vibv(u)=z(u)+vt(2*i+u) en 
!la subrutina mallaaguassomeras.
!Si no tiene sentido una CI tipo plano con altura constante (lo tiene si z(CCv)>CCH) será difícil darla. Para la posterior selección del dominio 
!se puede usar Manning y para el cálculo estacionario no hará falta definir nada si ap='si'. 
vibv(u)=z(u)-0.1_8 
if (z(u).le.33.48_8)then
vib(2*i+u)=34.48_8  	
else
vib(2*i+u)=z(u)+1.48_8
endif
!En ejemplos utilizados (planos con altura constante):
!Un plano de altura constante (tanto para altura como para nivel freático). 										 
!vib(2*i+u)=34.48_8
!vibv(u)=34.48_8
!Un plano de altura constante si z<altura y una superficie con altura mayor que ese valor en otro caso (tanto menor que la del terreno
!tanto mayor su cota). Visto que funciona mejor que un plano con altura constante si las CC de caudal subterráneo son nulas y se aplica lluvia.
!50.0_8 es un parámetro de ajuste de la superficie.
!if (z(u).lt.34.48_8)then	 
!vib(2*i+u)=34.48_8
!vibv(u)=34.48_8 	
!else
!vib(2*i+u)=z(u)-(z(u)-34.48_8)/50.0_8	 
!vibv(u)=z(u)-(z(u)-34.48_8)/50.0_8
!endif
enddo					

!Selección de subdominios superficiales bien desde la CI o la aproximación inicial (caso estacionario), bien desde el número de manning. 
!---------------------------------------------------------------------------------------------------------------------------------------  
!En caso de usar el modelo superficial con N-S 2D o el modelo subterráneo, no se realizará esta selección.
if ((modelo.eq.'conjunto').or.((modelo.eq.'superficial').and.(navier.eq.'no'))) then
do u=1,i
vb(u)=vib(u)
vb(i+u)=vib(i+u)
vit(u)=vib(u) 
vit(i+u)=vib(i+u)
vitv(u)=0.0_8
eval(u)=sqrt(3.0_8) 
!Selección de dominio inicial con la CI o aproximación inicial.
!if (vibv(u).lt.z(u)) then
!Selección de dominio inicial con el número de manning (selecc como superf zona con manning 0.05).	 !es
if (ma(u).ne.0.05_8) then  								                                             !es
evol(u)=0.0_8
else
evol(u)=sqrt(3.0_8)
endif
enddo
!Sobreescritura del coeficiente de Manning tras seleccionar el dominio con él (hasta enddo). !es
do u=1,i		 
ma(u)=0.0_8		
enddo			
!Obtención de una aproximación a partir de NS2D con selección del dominio posterior a su resolución	(bueno si velocidades nulas en aprox inicial o CI)
ap='si'

!Selección del subdominio 																		  
open(unit=4,file='C:\mallainicial.txt',status='old') 
read(4,27)a	 
 do u=1,j
 read(4,26)no(1),no(4),no(2),no(5),no(3),no(6)   
   if ((evol(no(1)).gt.0.0_8).and.(evol(no(2)).gt.0.0_8).and.(evol(no(3)).gt.0.0_8)) then
   do ui=1,3
   vit(2*i+no(ui))=vib(2*i+no(ui))-z(no(ui))
   vb(2*i+no(ui))=vib(2*i+no(ui))
   enddo
   else  
   do ui=1,3
   vitv(no(ui))=vibv(no(ui)) 
   eval(no(ui))=0.0_8
   eval(no(3+ui))=0.0_8
   enddo
   endif 
 enddo
close(4) 
!Se observa si existe flujo superficial en todo el dominio o si no existe flujo superficial en todo el dominio.
!Si se resuelve el modelo conjunto se resolverá (inicialmente) sólo la ecuación de aguas someras o sólo la ecuación subterránea.
if (minval(eval).eq.sqrt(3.0_8)) then 
write(6,*)'...Solo flujo superficial para esa condicion inicial o aproximacion'
elseif (maxval(eval).eq.0.0_8) then
write(6,*)'...No hay flujo superficial para esa condicion inicial o aproximacion'		
endif
endif

!Elección de modelización estacionaria (esquema implícito) o transitoria (esquema implícito o semi-implícito)
!------------------------------------------------------------------------------------------------------------
!Se construirán las cajas de masa para las ecuaciones de aguas someras y para la ecuación subterránea para el caso transitorio.
!También se impone condición de lluvia por pantalla, que es considerada para el modelo superficial (si aguas someras y formulación BG) y subterráneo.
if (imp.eq.'no')then
sino='si'
else
write(6,*)'Quieres obtener la solucion en tiempos en vez de'
write(6,*)'obtener la solucion estacionaria? (si/no)'
read(5,*)sino
endif
do while (nt.eq.0)
 if (sino.eq.'si') then	
 !Simulación transitoria
  write(6,*)' '
  write(6,*)'t=',Ata
  write(6,*)'Escribe el intervalo de tiempo total de simulacion en segundos:'
  read(5,*)tiempo
  write(6,*)'Escribe el incremento de tiempo de resolucion en segundos:'	   
  read(5,*)At
  !Se permite la impresión cada cierto número de incrementos de simulación para evitar generar archivos de mucho tamaño (por ejemplo, imprimir cada 
  !2 permitirá generar 5 archivos en vez de 10 al simular tiempo=10 días con At=1 día). Así, habrá incrementos de tiempo en que no se imprima.
  write(6,*)'Cada cuantos incrementos se imprime?'
  read(5,*)frec
  tac=tiempo+tac
  nt=idint((tac-Ata)/At)
  if (((navier.eq.'no').and.(modelo.ne.'subterraneo')).or.(modelo.eq.'subterraneo')) then
  !El volumen de agua por lluvia será multiplicado por At (para cada modelo).
  write(6,*)'Introduce un valor de intensidad media de lluvia'
  write(6,*)'durante el intervalo de simulacion en mm/d (habitual 0-100)'
  endif
 else 
 !Simulación estacionaria
  nt=1
  if (((navier.eq.'no').and.(modelo.ne.'subterraneo')).or.(modelo.eq.'subterraneo')) then
  !El volumen de agua por lluvia sería infinito para cualquier intensidad.  
  write(6,*)'Introduce un valor de intensidad de lluvia en mm/d'
  write(6,*)'Debe ser constante en el tiempo (habitual 0-100)'
  endif
 endif
 if	(((navier.eq.'no').and.(modelo.ne.'subterraneo')).or.(modelo.eq.'subterraneo')) then
 read(5,*) qll
 do u=1,i
 !Intensidad de lluvia en m/s (1mm=1L/m2)
 ql(u)=qll/(86400*1000)	
 enddo
 endif
do ut=1,nt
!Problema no lineal superficial:
!En aguas someras no es posible tener valores nulos o negativos de calado como valores de la iteración anterior. Si aparecen después la condición 
!seco-mojado limita el dominio de resolución.
!Para transitorio se utiliza la solución en instante anterior en cada paso de tiempo para comenzar Picard o Newton.
!Se utiliza la condición inicial en el primer paso si aguas someras en modelo superficial o conjunto (sino la trivial). 
!Para estacionario se utiliza la aproximación inicial si aguas someras en modelo superficial o conjunto (sino la trivial).
do u=1,3*i
vt(u)=vit(u)     
enddo
!Problema no lineal subterráneo:
!En ec. subterránea tampoco es posible tener valores nulos o negativos de calado como valores de la iteración anterior. Si aparecen después, la condición 
!de espesor mínimo permite calcular coeficientes que no dan problema.
!Para transitorio se utiliza la solución en instante anterior en cada paso de tiempo para comenzar Picard.
!Se utiliza la condición inicial en el primer paso en modelo conjunto (sino la solución 0.1).
!Para estacionario se utiliza la aproximación inicial en modelo conjunto (sino la solución 0.1).
do u=1,i
vtv(u)=vitv(u)    
enddo
													   	  
!Gestión de las ecuaciones de N-S 2D, aguas someras y agua subterránea.
!----------------------------------------------------------------------
!Resolución de los modelos superficial, subterráneo y conjunto.  
!Respecto a 'vn' y 'vnv', en caso de tratarse del modelo conjunto (vn,vnv) o del modelo superficial (vn) y existe contorno móvil (aplicación de condición 
!seco-mojado o de la condición similar para el modelo conjunto), pueden ser diferentes en cada iteración no-lineal de aguas someras. Cada uno llevará los 
!nodos cuadráticos que tienen la malla superficial y la subterránea respectivamente.
!En las ecuaciones de aguas someras con Manning=0 y Re1 se obtendrán soluciones similares a si se usa Manning=n y Re2 con Re2>Re1.

!Modelo superficial.
if (modelo.eq.'superficial') then
 tol=1e-5
 it=0
 !La condición seco mojado se aplica dentro de la subrutina aguassomeras: a partir de la segunda iteración 
 !no-lineal se utiliza el resultado de flujo superficial de la iteración anterior con velocidades nulas en el contorno móvil.
 !Permite añadir elementos al dominio superficial o desechar elementos del dominio superficial. 
 if ((navier.eq.'no').and.(itt.eq.0)) then 
 !Necesario usar la subrutina mallaaguassomeras si se selecciona dominio previamente con la CI o la aproximación inicial (aquí también se selecciona
 !la malla para flujo superficial, algo que después hace la condición seco-mojado) 
 call mallaaguassomeras (i,j,z,vtv,velx,vely,eval,vo,vt)  
  !Aproximación inicial con las ecuaciones N-S 2D. Necesario haber usado un esquema implícito ya que que N-S 2D sólo está programado para implícito. 
  if (ap.eq.'si') then
  write(6,*)'Ojo, primera/as it con NS2D'                                                                                              
  !Podría usarse un dominio menor al de toda la malla al haberse hecho una selección previa (buscar un Re que permita calados en el dominio 
  !seleccionado.
  nu=0.5_8 
  !Seleccionar el número de iteraciones (en el código). 
  !Estas iteraciones se hacen con Picard ya que N-S 2D sólo está programado para Picard.
  !Se hará una selección tras cada resolución con la condición seco-mojado.  
  !do u=1,6                                                                                                                       
  call aguassomeras(modelo,'si',ap,tol,i,j,x,y,z,it,itt,nu,ma,sino,bcg,ten,'no',At,vn,vo,eval,ql,vib,vt,vb,est,imp,vit)           
  !enddo                                                                                                                          
  nu=mu
  endif
 endif 
 !En aguas someras se aplica el método de Picard y el de Newton y se hará una selección tras cada resolución con la condición seco-mojado.
 !Esquema semi-implícito - iteración para un incremento de tiempo.
 if (imp.eq.'no')then
 call aguassomeras(modelo,navier,ap,tol,i,j,x,y,z,it,itt,nu,ma,sino,bcg,ten,newton,At,vn,vo,eval,ql,vib,vt,vb,est,imp,vit)
 write(6,*)'t=',Ata+At
 else
 !Esquema implícito - inicio de iteraciones para un incremento de tiempo o para el problema estacionario.
 do while (tol.eq.1e-5)
 call aguassomeras(modelo,navier,ap,tol,i,j,x,y,z,it,itt,nu,ma,sino,bcg,ten,newton,At,vn,vo,eval,ql,vib,vt,vb,est,imp,vit)  
 enddo
 write(6,*)'Convergencia modelo aguassomeras/navierstokes en iteracion:',it
  if (sino.eq.'si')then
  write(6,*)'At=',Ata+At
  endif
 endif
 !En ejemplos anteriores:
 !Resultados velocidades en linea central del cuadrado para Cavity flow
 !open(unit=11,file='C:\soluciones linea central.txt',status='unknown')
 !228 format (2x,F10.5,2x,F10.5)
 !do u=1,i
 !if (x(u).eq.200) then
 !write(11,228) vt(u),y(u)/400
 !endif
 !enddo
 !close(11)
 call vorticidad (i,j,x,y,nu,vt,te,vor)

!Modelo subterráneo.
elseif (modelo.eq.'subterraneo') then 
 !Se usa toda la malla.
 !No se aplica ninguna condición que permita seleccionar elementos o desechar elementos del dominio subterráneo.
 !Se aplica el método de Picard.
 !Esquema implícito - Dentro de la subrutina aguassubterránea se hacen iteraciones para un incremento de tiempo o para el problema estacionario.
 call aguassubterranea(i,j,x,y,zp,nd,kix,kiy,ag,sino,bcg,At,vn,vov,eval,ql,qb,qx,qy,vibv,vtv,velx,vely) 
 if (sino.eq.'si')then
 write(6,*)'At=',Ata+At
 endif

!Modelo conjunto. 
!Se decide si existe flujo superficial y subterráneo, y si existen ambos se aplican alternativamente las ecuaciones de 
!aguas someras y la ecuación de flujo subterráneo. En este caso se realiza una iteración de las primeras y se resuelve la segunda (en cada 
!iteración del modelo conjunto). Se realizan estas iteraciones hasta convergencia del modelo superficial, y se hace para 
!cada incremento de tiempo (partiendo de la misma condición inicial) si se busca una solución transitoria.
!Se resuelve primero el superficial aunque no haya condiciones de contorno de flujo superficial lo que puede ocurrir si se ha seleccionado previamente
!con la CI-aprox una zona de agua embalsada (en cuyo caso ya no se habrán definido condiciones de este tipo en el contorno de todo el dominio, y como
!no se están usando CC variables en el tiempo siempre se usarán siempre éstas en caso de simular un transitorio...).
elseif (modelo.eq.'conjunto') then   
 tol=1e-6
 it=0
 !La condición similar a la condición seco mojado se aplica dentro de la subrutina aguassomeras: a partir de la segunda iteración 
 !no-lineal se utiliza el resultado de flujo superficial de la iteración anterior con velocidades subterráneas en el contorno móvil.
 !Permite añadir elementos al dominio superficial o desechar elementos del dominio superficial (y generar flujo subterráneo aislado).
 !No se aplica ninguna condición en la ecuación subterránea que permita seleccionar elementos o desechar elementos del dominio 
 !subterráneo (ni generar flujo superficial aislado), modificándose éste con la modificación del dominio superficial.
 !Inicio de iteraciones conjuntas: 
  do while (tol.eq.1e-6)
  if (maxval(eval).eq.0.0) then 
  !Toda la solución previa es subterránea (de acuerdo a la selección previa ó a la última solución), por lo que se resuelve implícitamente el 
  !subterráneo hasta convergencia (Picard) y se tiene la solución conjunta para este incremento de tiempo, o la solución conjunta estacionaria.	       
   call aguassubterranea(i,j,x,y,zzp,nd,kix,kiy,ag,sino,bcg,At,vn,vov,eval,ql,qb,qx,qy,vibv,vtv,velx,vely)
   !Las siguientes dos líneas son equivalentes a escribir: exit
   tol=1.0
   cycle   
  else
  !La solución previa es toda superficial o superficial y subterránea. 	 
  !Siempre necesario usar la subrutina mallaaguassomeras previamente para obtener las velocidades subterráneas a aplicar en el contorno móvil (aquí
  !también se selecciona la malla para flujo superficial algo que ahora no hace la condición similar a la condición seco-mojado).
  call mallaaguassomeras (i,j,z,vtv,velx,vely,eval,vo,vt)
   if (minval(eval).eq.0.0) then
   !Hay solución superficial en al menos un elemento (si sol. superficial sólo en un nodos aislados no se calcula con la ec. aguas someras).
   !Aproximación inicial con las ecuaciones N-S 2D (las ecuaciones de aguas someras también se resolverán implícitamente).
   if ((it.eq.0).and.(itt.eq.0).and.(ap.eq.'si'))then
   write(6,*)'Ojo, primera/as it con NS2D'
   !Podría usarse un dominio menor al de toda la malla al haberse hecho una selección previa (buscar un Re que permita calados en el dominio 
   !seleccionado.                                                                                                  
   nu=0.5_8	   
   !Seleccionar el número de iteraciones (en el código). 
   !Estas iteraciones se hacen con Picard ya que N-S 2D sólo está programado para Picard.
   !Se hará una selección tras cada resolución con la condición seco-mojado (no la similar a la condición seco-mojado).
   !do u=1,6                                                                                                                       
   call aguassomeras('superficial ','si',ap,tol,i,j,x,y,z,it,itt,nu,ma,sino,bcg,ten,'no',At,vn,vo,eval,ql,vib,vt,vb,est,imp,vit)   
   !enddo                                                                                                                          
   nu=mu
   endif
   !Se aplica implícitamente una iteración de la ecuación de aguas someras (Picard o Newton) para este incremento de tiempo, o para el cáculo de la 
   !solución conjunta estacionaria.	Se hará una selección tras la resolución con la condición similar a la condición seco-mojado.
   call aguassomeras(modelo,navier,ap,tol,i,j,x,y,z,it,itt,nu,ma,sino,bcg,ten,newton,At,vn,vo,eval,ql,vib,vt,vb,est,imp,vit) 	      
   endif
  endif
  if (maxval(eval).eq.0.0) then
  !Toda la solución previa es superficial. 
  !No se calcula solución subterránea y se cambia eval para que no indique posteriormente que toda el agua es subterránea en la siguiente
  !iteración conjunta.
   do u=1,i
   eval(u)=sqrt(3.0_8)
   enddo
  else
  !La solucion previa en subterránea o superficial y subterránea
  !Siempre necesario usar la subrutina mallasubterranea previamente para obtener las alturas superficiales a aplicar en el contorno móvil 
  !(que irán en vov(i)). Aquí también se selecciona la malla para flujo subterráneo, definida por el dominio superficial (la que malla que queda).   						 
   call mallasubterranea (i,j,z,zp,zzp,vt,eval,qx,qy,vov)
   if (minval(eval).eq.0.0) then
   !Hay solución subterránea en al menos un elemento (si sol. subterránea sólo en nodos aislados no se calcula con la ec. subterránea).
   call aguassubterranea(i,j,x,y,zzp,nd,kix,kiy,ag,sino,bcg,At,vn,vov,eval,ql,qb,qx,qy,vibv,vtv,velx,vely) 
   endif
  endif	   
  enddo
  write(6,*)'Convergencia modelo conjunto en iteracion:', it
  if (sino.eq.'si')then
  write(6,*)'At=',Ata+At
  endif	 
endif 

!Arreglo de la solución (interpolación en nodos cuadráticos y cálculo de variables conjuntas para el modelo conjunto).
!---------------------------------------------------------------------------------------------------------------------
if ((modelo.eq.'superficial').or.(modelo.eq.'conjunto')) then	  
!Para la próxima iteración no-lineal ec. superficial en el siguiente paso de tiempo (transitorio) para el esquema implícito (ambos modelos) o para 
!la próxima iteración no-lineal en el siguiente paso de tiempo para el esquema semi-implícito (modelo superficial), se considera la solución 
!superficial obtenida (necesarios calados en ella):
do u=1,3*i
vit(u)=vt(u)
enddo
!Las alturas en los nodos esquina calculadas debajo difieren de las calculadas en 'vb' en que llevan las cotas del terreno en el subdominio 
!subterráneo (modelo conjunto). Así, saldrá la cota de la lámina de agua si hay calado y la cota del terreno si el calado es nulo.  
do u=1,i
vt(2*i+u)=vt(2*i+u)+z(u)
vx(u)=0.0_8
vy(u)=0.0_8
ht(u)=0.0_8
enddo
open(unit=1,file='C:\malla.txt',status='old')
!Operaciones hechas en los elementos del dominio superficial.
!Las alturas irán en la variable conjunta ht (se completará más adelante en el caso del modelo conjunto).
!Las velocidades irán en las variables conjuntas vx,vy (se completará más adelante en el caso del modelo conjunto). 
read(1,21)j
!Necesario leer j con los elementos superficiales (j será el número de elementos subterráneos si se ha usado el modelo conjunto).
read(1,'(A)')a
do u=1,j
 read(1,26)no(1),no(4),no(2),no(5),no(3),no(6) 
   do ui=1,3
   b=ui+1-sb(ui)
   c=ui+2-sb(ui+1)   
   !Interpolacion lineal postproceso de alturas y Manning en los nodos cuadráticos.
   vt(2*i+no(ui+3))=(vt(2*i+no(ui))+vt(2*i+no(b)))/2.0_8  
   ma(no(ui+3))=(ma(no(ui))+ma(no(b)))/2.0_8
   !Se escriben los valores de las variables para flujo superficial en las variables conjuntas.
   vx(no(ui))=vt(no(ui))	 
   vx(no(ui+3))=vt(no(ui+3))
   vy(no(ui))=vt(i+no(ui))
   vy(no(ui+3))=vt(i+no(ui+3))
   ht(no(ui))=vt(2*i+no(ui))
   ht(no(ui+3))=vt(2*i+no(ui+3))
   enddo
enddo
close(1)
endif
if ((modelo.eq.'subterraneo').or.(modelo.eq.'conjunto')) then
!Para la próxima iteración no-lineal ec. subterránea en el siguiente paso de tiempo (transitorio) para el esquema implícito (ambos modelos), se 
!considera la solución subterránea obtenida (necesarias alturas en ella):
do u=1,i
vitv(u)=vtv(u)
enddo
open(unit=3,file='C:\mallasub.txt',status='old')
!Operaciones hechas en los elementos del dominio subterráneo.
read(3,21)j
read(3,'(A)')a
do u=1,j
read(3,26)no(1),no(4),no(2),no(5),no(3),no(6) 
 do ui=1,3
 b=ui+1-sb(ui)
 c=ui+2-sb(ui+1)
 !Interpolación lineal postproceso de la velocidad subterránea, la altura (nivel freático) y las propiedades en los nodos cuadráticos.
 vtv(no(ui+3))=(vtv(no(ui))+vtv(no(b)))/2.0_8
 kix(no(ui+3))=(kix(no(ui))+kix(no(b)))/2.0_8
 kiy(no(ui+3))=(kiy(no(ui))+kiy(no(b)))/2.0_8
 ag(no(ui+3))=(ag(no(ui))+ag(no(b)))/2.0_8
 nd(no(ui+3))=(nd(no(ui))+nd(no(b)))/2.0_8
 velx(no(ui+3))=(velx(no(ui))+velx(no(b)))/2.0_8
 vely(no(ui+3))=(vely(no(ui))+vely(no(b)))/2.0_8
 !Se completan las variables conjuntas escribiendo ahora en ellas los valores de las variables para flujo subterráneo.
 vx(no(ui))=velx(no(ui))	 
 vx(no(ui+3))=velx(no(ui+3))
 vy(no(ui))=vely(no(ui))
 vy(no(ui+3))=vely(no(ui+3))
 ht(no(ui))=vtv(no(ui))
 ht(no(ui+3))=vtv(no(ui+3))
  !Se darán valores de calado en ciertos nodos cuadráticos que no tienen ecuación superficial     !cñ  
  !Para las condiciones puestas, no(ui) pertenecerá al dominio superficial, no(b) al subterráneo y no(ui+3) será el nódo cuadrático entre ellos.   
  !Con valores de nivel freático:
  !if ((vtv(no(ui+3)).gt.z(no(ui+3))).and.(vt(2*i+no(ui)).gt.z(no(ui))).and.(vt(2*i+no(b)).eq.z(no(b)))) then
  !vt(2*i+no(ui+3))=vtv(no(ui+3))
  !elseif ((vtv(no(ui+3)).gt.z(no(ui+3))).and.(vt(2*i+no(ui)).eq.z(no(ui))).and.(vt(2*i+no(b)).gt.z(no(b)))) then
  !vt(2*i+no(ui+3))=vtv(no(ui+3))
  !endif
  !Con valores de altura de nodos cercanos del dominio superficial (altura supuesta constante en el flujo superficial):
  if ((vt(2*i+no(ui)).gt.z(no(ui+3))).and.(vt(2*i+no(ui)).gt.z(no(ui))).and.(vt(2*i+no(b)).eq.z(no(b)))) then
  vt(2*i+no(ui+3))=vt(2*i+no(ui))
  !Para que la variable conjunta de la altura (ht) coincida con la variable altura superficial o cota del terreno (vt)
  !(con afección en caso de aplicar el modelo conjunto): ht(no(ui+3))=vt(2*i+no(ui+3)) 
  elseif ((vt(2*i+no(b)).gt.z(no(ui+3))).and.(vt(2*i+no(ui)).eq.z(no(ui))).and.(vt(2*i+no(b)).gt.z(no(b)))) then
  vt(2*i+no(ui+3))=vt(2*i+no(b))
  !Para que la variable conjunta de la altura (ht) coincida con la variable altura superficial o cota del terreno (vt): ht(no(ui+3))=vt(2*i+no(ui+3))
  endif
 enddo
enddo
close(3)
endif
!Otras variables que se pueden calcular postproceso:  !ci
!Velocidad vertical para la ecuación de aguas someras a una zh (siguiendo a Raquel Taboada, capítulo 5).
!Se utilizaría la subrutina pendientes y mod=1 para tener pe(u)=dz/dx y pe(u+i)=dz/dy. Habría que calcular du/dx y dv/dy.
!mod=1.0_8
!call pendientes (i,j,x,y,z,mod,pe)
!do u=1,i
!w(zh)=vt(u)*pe(u)+vt(i+u)*pe(i+u)-(zh-z(u))*(du/dx+dv/dy)
!enddo
!Módulo de la velocidad lineal para la ecuación para flujo subterráneo. Es la vel. real (vel. subterránea/porosidad efectiva), diferente a la 
!vel. subterránea o vel. de Darcy (caudal/secc. total) calculada, y es mayor ya que la sección es menor al haber poros. 
!En el modelo conjunto se usa la velocidad Darcy para generar el caudal que tiene lugar (al distribuirla en el contorno). No hay zonas del contorno
!sin entrada de flujo debido al medio poroso. 
!do u=1,i
!vell(u)=sqrt(velx(u)**2.0_8+vely(u)**2.0_8)/nd(u)
!enddo

!Escritura de la solución en fichero:
!------------------------------------
!Se escribe en el fichero solucionaguassomeras-subterraneo con formatos preparados para el programa Tecplot.	  
write(6,*)' '
if (sino.eq.'si') then
!Escritura para el esquema implícito con caso transitorio o para el esquema	semi-implícito (permite posterior animación).
Ata=Ata+At
itt=itt+1
 !Se decide si se imprime la solución (para algunos incrementos de tiempo podría no hacerse).
 if (ut.eq.frec*ni) then
 !Se pasan números a formato caracter y se concatenan cadenas de caracteres. 
 write(Atee,'(i10)')idint(Ata)
 Atai='"solucion para t='//Atee
 Atai=Atai(1:27)//'s"'
 endif					 
else
!Escritura de una sola solución para el esquema implícito con caso estacionario.
Atai='"solucion estacionaria unica"'
endif
if ((ut.eq.frec*ni).or.(sino.eq.'no')) then	 	  
write(8,22)Atai
endif					 

!Escritura de las variables representativas (con doble precisión) para cada modelo.
if (modelo.eq.'superficial') then
 !Modelo superficial. Para el siguiente paso de tiempo ec. superficial (transitorio) para el esquema implícito o para el siguiente paso de tiempo 
 !para el esquema semi-implícito, se considera como solución en el instante anterior la obtenida (necesarias alturas en ella):
 do u=1,3*i
 vib(u)=vt(u) 
 enddo
 !Se decide si se imprime la solución (para algunos incrementos de tiempo podría no hacerse).
 if ((ut.eq.frec*ni).or.(sino.eq.'no')) then	 
 open(unit=4,file='C:\mallainicial.txt',status='old')
 !Impresión de valores de coordenadas x e y en metros, de altura de la lámina superficial-terreno en metros, de calado (incluido en NS-2D) en metros, 
 !de velocidad superficial (en direcciones x e y) en m/s, del módulo de la velocidad, de tensión (en direcciones x e y) en kilopascales, 
 !de tensión tangencial xy en kilopascales (kPa=1000Pa), y de vorticidad en 1/s.
 read(4,21)j
 write(8,23)i,4*j
 write(8,'(A)')'DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE )'
 do u=1,i
 te(u)=te(u)/1000.0_8
 te(i+u)=te(i+u)/1000.0_8
 te(2*i+u)=te(2*i+u)/1000.0_8
 write(8,38)x(u),y(u),vt(2*i+u),vt(2*i+u)-z(u),vt(u),vt(i+u),sqrt(vt(u)**2+vt(i+u)**2),te(u),te(i+u),te(2*i+u),vor(u)
 enddo
 endif					 
elseif (modelo.eq.'subterraneo') then
 !Modelo subterráneo. Para el siguiente paso de tiempo ec. subterránea (transitorio) se considera como solución en el instante anterior la obtenida
 !(necesarias alturas en ella):
 do u=1,i
 vibv(u)=vtv(u) 
 enddo
 !Se decide si se imprime la solución (para algunos incrementos de tiempo podría no hacerse).
 if ((ut.eq.frec*ni).or.(sino.eq.'no')) then	 
 open(unit=4,file='C:\mallasubinicial.txt',status='old')
 !Impresión de valores de coordenadas x e y en metros, de nivel freático en metros, de espesor freático en metros, de velocidad de Darcy 
 !(en direcciones x e y) en m/s, y del módulo de la velocidad.
 read(4,21)j
 write(8,23)i,4*j
 write(8,'(A)')'DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE )'
 do u=1,i
 write(8,39)x(u),y(u),vtv(u),vtv(u)-zp(u),velx(u),vely(u),sqrt(velx(u)**2+vely(u)**2)
 enddo
 endif	                 
else
 !Modelo conjunto. Para el siguiente paso de tiempo ec. superficial y subterránea (transitorio) se considera como solución en el instante anterior la 
 !solución conjunta obtenida (variables conjuntas). También es necesaria una solución con las alturas en ella:
 do u=1,i
 vib(u)=vx(u)
 vib(i+u)=vy(u)
 vib(2*i+u)=ht(u)
 vibv(u)=ht(u)
 enddo
 if ((ut.eq.frec*ni).or.(sino.eq.'no')) then
 !Se decide si se imprime la solución (para algunos incrementos de tiempo podría no hacerse).	 
 open(unit=4,file='C:\mallainicial.txt',status='old')
 !Sólo aquí se utilizan de las variables conjuntas ht, vx y vy.
 !Impresión de valores de coordenadas x e y en metros, de altura de la lámina superficial-nivel freático (variable conjunta) en metros, de altura de la 
 !lámina superficial-terreno en metros, de calado (incluido en NS-2D) en metros, de espesor freático en metros, de velocidad superficial-velocidad de 
 !Darcy (en direcciones x e y, variables conjunta) en m/s, y del módulo de la velocidad.
 read(4,21)j	
 write(8,23)i,4*j
 write(8,'(A)')'DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE )'
 do u=1,i
 !Hay que tener en cuenta que el modelo conjunto no calcula flujo subterráneo bajo el superifical.
 !Por defecto, siempre habrá calado nulo donde no se calcula flujo superficial (vt-z=vt original). Pero no ocurrirá lo mismo con el espesor freático
 !donde no se calcula flujo subterráneo. Aunque podría darse z(u)-zp(u), para ser consecuente se dará espesor freático nulo.
 !Ello permitirá reconocer hasta donde resuelve cada sub-modelo al representar estas variables.
 if (vtv(u).ne.0.0_8)then
 !En caso de calcularse flujo subterráneo se imprimen todas las variables (el calado será nulo).
 write(8,40)x(u),y(u),ht(u),vt(2*i+u),vt(2*i+u)-z(u),vtv(u)-zp(u),vx(u),vy(u),sqrt(vx(u)**2+vy(u)**2)
 else
 !En caso de no calcularse flujo subterráneo se da espesor freático nulo. 
 write(8,40)x(u),y(u),ht(u),vt(2*i+u),vt(2*i+u)-z(u),0.0_8,vx(u),vy(u),sqrt(vx(u)**2+vy(u)**2)
 endif
 enddo		
 endif  	                 
endif

!Escritura de la malla.
if ((ut.eq.frec*ni).or.(sino.eq.'no')) then	 
 ni=ni+1					 
 write(8,'(A)')' '
 read(4,'(A)')a
 do u=1,j
 read(4,26)no(1),no(4),no(2),no(5),no(3),no(6)
 !Para tecplot cada elemento ha de ser lineal por lo que cada elemento supone cuatro a representar. En total 4*j elementos.
 !Por tanto se apreciará una malla equivalente al representarla.
 write(8,41) no(1),no(4),no(6)
 write(8,41) no(4),no(2),no(5)
 write(8,41) no(5),no(3),no(6)
 write(8,41) no(4),no(5),no(6)
 enddo
 close(4)
 write(6,*)'Fichero solucion aguassomeras-subterraneo creado'
endif                   
enddo

!Se decide si se termina la modelización:
!----------------------------------------
if (sino.eq.'si') then
write(6,*)'Quieres obtener la simulacion para un tiempo mayor? (si/no)'
read(5,*)sino
  if (sino.eq.'si') then
  nt=0
  ni=1				    
  endif
endif
enddo

!Escritura de propiedades (constantes para todos los incrementos de tiempo):
!---------------------------------------------------------------------------
!Se evita escribir datos innecesario en caso de utilizar el caso transitorio con el esquema implícito o el esquema semi-implícito. 
!Se escribe en el fichero solucionpropiedades con formatos preparados para el programa Tecplot.
if (modelo.eq.'superficial') then
 open(unit=4,file='C:\mallainicial.txt',status='old')
 read(4,21)j
 !Valores de coordenadas en metros (con cota del terreno), y de Manning,
 write(10,'(A)')'TITLE = "Propiedades para aguassomeras"'
 write(10,'(A)')'VARIABLES =  X, Y, Z, n'
 write(10,24)i,4*j
 write(10,'(A)')'DT=(DOUBLE DOUBLE DOUBLE DOUBLE )'
 do u=1,i
 write(10,48)x(u),y(u),z(u),ma(u)
 enddo 
elseif (modelo.eq.'subterraneo') then
 open(unit=4,file='C:\mallasubinicial.txt',status='old')
 read(4,21)j
 write(10,'(A)')'TITLE = "Propiedades para flujo subterraneo"'
 !Valores de coordenadas en metros (con cota del sustrato), de conductividad (en direcciones x e y) en m/s, de ángulo de anisotropía en grados, 
 !y de porosidad (adimensional).
 write(10,'(A)')'VARIABLES =  X, Y, ZP, Kix, Kiy, Angº, nd'
 write(10,24)i,4*j
 write(10,'(A)')'DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE )'
 do u=1,i
 write(10,49)x(u),y(u),zp(u),kix(u),kiy(u),ag(u)*180.0_8/3.14159_8,nd(u)
 enddo
else
 open(unit=4,file='C:\mallainicial.txt',status='old')
 read(4,21)j
 write(10,'(A)')'TITLE = "Propiedades para aguassomeras y subterraneo"'
 !Valores de coordenadas en metros (con cota del terreno y cota del sustrato), de conductividad (en direcciones x e y) en m/s, de ángulo 
 !de anisotropía en grados sexagesimales, de porosidad (adimensional) y de Manning.
 write(10,'(A)')'VARIABLES =  X, Y, Z, ZP, Kix, Kiy, Angº, nd, n'
 write(10,24)i,4*j
 write(10,'(A)')'DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE )'
 do u=1,i
 write(10,50)x(u),y(u),z(u),zp(u),kix(u),kiy(u),ag(u)*180.0_8/3.14159_8,nd(u),ma(u)
 enddo
endif
!Escritura de la malla
write(10,'(A)')' '
read(4,'(A)')a
do u=1,j
read(4,26)no(1),no(4),no(2),no(5),no(3),no(6)
write(10,41) no(1),no(4),no(6) 
write(10,41) no(4),no(2),no(5) 
write(10,41) no(5),no(3),no(6) 
write(10,41) no(4),no(5),no(6) 
enddo
close(4)

!Eliminación de los ficheros generados:
!--------------------------------------
if ((modelo.eq.'superficial').or.(modelo.eq.'conjunto')) then
open(unit=1,file='C:\malla.txt',status='old')
close(1,status = 'delete')
endif
if ((modelo.eq.'subterraneo').or.(modelo.eq.'conjunto')) then
open(unit=3,file='C:\mallasub.txt',status='old')
close(3,status = 'delete')
endif

deallocate(x,y,z,zp,zzp,vt,vit,vib,vb,vo,eval,vov,vtv,vitv,vibv,kix,kiy,ag,nd,ma,qx,qy,qb,ql,velx,vely,vx,vy,ht,evol,te,vor)
end

!----------------------------------------------------------------------------------------------------------------------------------------------------------
!Subrutinas de gestión para el cálculo de la solución de las ecuaciones.
!----------------------------------------------------------------------------------------------------------------------------------------------------------
!Subrutina AGUASSOMERAS (tres ecuaciones, 2 dinámicas y 1 de continuidad). En esta subrutina se calculan las velocidades, la altura de la lámina de agua. 
!Se construye el sistema con la matriz de rigidez (de cajas en las que irán A, B,...). Se calculan las alturas (respecto a z=0) en nodos esquina y las 
!velocidades en todos los nodos. Se calculan los calados postproceso. 
!Se pueden resolver por el método de Picard la ecuación de aguas someras o la ecuación de N-S 2D.
!Se puede utilizar el método de Newton para resolver la ecuación de aguas someras, en cuyo caso se aplicarán 5 iteraciones previas de Picard.
!Sólo se hará una iteración no-lineal y se construirán todas las cajas cada vez (la posible afección de la condición seco-mojado obliga 
!a ello aunque las cajas no contengan coeficientes no lineales).  
!El vector inicial (it=0, itt=0) para resolver la no linealidad estará en vt y vb (aprox inicial en estacionario, CI u otro valor en transitorio) y tendrá
!los valores Dirichlet prescritos (variables del sistema) de la altura de la lámina de agua y de velocidad para el cálculo de las integrales de contorno.
!Al realizarse la primera iteración de Picard (it=0, itt=0) se pueden resolver algunas iteraciones de las ecuaciones de N-S 2D al resolver 
!la ecuación de aguas someras. Así se calcula un buen vector inicial para resolver en la siguiente iteración las ecuaciones de aguas someras.
!Siempre habrá un calado, altura y velocidad nulos donde no se calcula flujo superficial (en vt y vb).
!'tol' gestiona la parada del bucle cuando dos soluciones consecutivas son suficientemente similares,
!'vn,ve,ven,vo,vu' son contadores de nodos cuadráticos, lineales o con condición de contorno,
!'veac,veoc,veuc' son vectores que contienen las integrales de contorno,'vef' es el término fuente de caudal por lluvia para flujo superficial,
!'fx,fy' son los vectores donde se calculan las integrales de la pendiente de fricción (cuantifican la tensión de fondo), 
!'fz' es un vector del mismo tipo que se genera al aplicar el método de estabilización,
!'sa' es un vector donde se almacenan los coeficientes no nulos de la matriz del sistema (formato MSR), 
!'ita' es el vector puntero con la posición de estos coeficientes (formato MSR),
!'cia,ca' son inicialmente análogos a 'ita,sa' al ser usados para guardar por separado las matrices de masa. Si se resuelve con 
!precondicionador diagonal (subrutina gradientesbiconjugados para formato MSR) 'ita,sa' llevarán la matriz del sistema a resolver y 'cia,ca' son desechados.
!En otro caso se dimensiona 'cja', se copiarán 'ita,sa' a 'ca,cia,cja' (subrutina dslubc para formato CSC) y 'ita,sa son desechados. Así,
!'ca' será el vector donde se almacenan los coeficientes no nulos de la matriz del sistema (formato CSC), 
!'cia,cja' serán los vectores punteros con la posición de estos coeficientes (formato CSC),
!'ndim' es el nº de coeficientes no nulos que hay en la matriz del sistema almacenados en los vectores 'ita,sa' (varia al formar el sistema a resolver),
!'inc' es la dimensión de la matriz del sistema cuadrada si se usasen elementos cuadráticos para alturas/calados y velocidades (nºecuaciones*i), 
!'c' es la dimensión de la matriz del sistema cuadrada para resolver el sistema (en cada iteración) tras reducir el orden con las CC y considerar los 
!elementos lineales para alturas/calados, 
!'ru,rb' son los vectores del residual y de la solución para el método de Newton,
!'vec,voc' forman el término independiente del sistema, 'vv,vt,vb' son los vectores donde se guarda la solución del sistema (en cada iteración), 
!'yn' permite el cálculo de una solución de aguas someras con calados similares sea como sea la cota del terreno,
!'nonzero' permite utilizar la solución anterior para aplicar el método de los gradientes biconjungados,
!'vth' lleva la solución del calado para la aplicación de la condición seco-mojado (subrutina nuevamalla).                                                                               
!-----------------------------------------------------------------------------------------------------------------------------------------------------------
subroutine aguassomeras(modelo,navier,ap,tol,i,j,x,y,z,it,itt,nu,ma,sino,bcg,ten,newton,At,vnin,vo,eval,ql,vib,vt,vb,est,imp,vit)
use allocatacion
integer*4, dimension(:),allocatable::vn
integer*4 i,j,u,uu,c,vnin(i),it,itt,inc,ndim,k			 
real*8, dimension(:),allocatable::vbdin,vecdin,vv,vec,voc,vu,vef,fx,fy,fz,vth,veac,veoc,veuc,ru,rb		  	 
real*8 x(i),y(i),z(i),ma(i),vt(3*i),vb(3*i),vit(3*i),vib(3*i),tol,nu,vo(3*i),eval(i),At,ql(i),del									  		
logical nonzero 
character modelo*12,sino*2,newton*2,bcg*2,ten*2,navier*2,yn*2,est*2,imp*2,ap*2
!Aquí sí hace un dimensionamiento dinámico propiamente dicho con 'vbdin,vecdin' (sería obligatorio definirlas así si fuesen variables globales). 
!Otro ejemplo son las variables globales ita,sa (usadas para resolver con precondicionador diagonal) ó cia,ca,cja (usadas para resolver con 
!precondicionador LU) cuya memoria se destruirá (según el caso) en tiempo de ejecución de esta subrutina (tras resolver con ellos).
allocate(vn(i),vv(3*i),vec(3*i),voc(3*i),vu(3*i),vef(i),fx(i),fy(i),fz(i),vth(i),veac(i),veoc(i),veuc(i),ru(3*i),rb(3*i))

30  format(2(3x,A1,'(',I5,')=',E15.8E2))	  

!Condición 'falsa' de calado normal. !cñ
!---------------------------------------
!Si yn='si' se genera una solución con calado aproximadamente idéntico en todos los nodos independiente de la cota del terreno.
!Sólo tenida en cuenta para newton='no' y est='no' y con sentido si navier='no' (aguas someras).
!Si yn='no' se resuelven las ecuaciones apropiadas.
yn='no'  
!Para la aplicación del esquema implícito o semi-implícito (el semi-implícito será de segundo orden y conllevará a menor error).
if (imp.eq.'no') then
del=2.0_8
else
del=1.0_8
endif

!Selección de en qué iteraciones se utilizará la solución anterior para el método de los gradientes biconjungados
!(aunque en caso de aplicar el método el Newton esto no será efectivo)
!----------------------------------------------------------------------------------------------------------------
if (it.eq.0)then
nonzero=.false.	
else
nonzero=.true.
endif

!Inizialización de variables:
!----------------------------
inc=3*i
do u=1,3*i
vv(u)=vt(u)
ru(u)=0.0_8
rb(u)=vb(u)
voc(u)=0.0_8
vu(u)=vo(u)
enddo
do u=1,i
fx(u)=0.0_8
fy(u)=0.0_8
fz(u)=0.0_8
vef(u)=0.0_8
veac(u)=0.0_8
veoc(u)=0.0_8
veuc(u)=0.0_8
enddo
!El vn original es ahora vnin. El nuevo vn será diferente si el dominio es menor al de toda la malla.
!En la ecuación de continuidad se eliminarán las ecuaciones relativas a todos los nodos cuadráticos y no sólo los del dominio superficial.
!Con vn se eliminarán los de este dominio y el resto se eliminarán con las CC al no considerar el resto de la malla.
do u=1,i
 if ((vnin(u).eq.1).and.(eval(u).eq.0.0)) then
 !Nodos cuadráticos y dominio superficial
 vn(u)=1
 else
 vn(u)=0
 endif
enddo
!Se dimensionan ita y sa (también cia y ca), de acuerdo al número total de coeficientes que se generarán con las matrices elementales.
!Además ita ya lleva referenciados el número de coeficientes que habrá por fila (en sus primeras 3*i+1 componentes) antes de dar el formato MSR.
!Cuando se calcule una matriz se almacenarán sus coeficientes directamente en estos vectores.
call dimvectas (i,j,newton,sino,vnin,est)
ndim=ita(inc+1)-1

!En caso de resolver de forma transitoria con incrementos de tiempo (esquema implícito considerando las matrices de masa o esquema semi-implícito).
!--------------------------------------------------------------------------------------------------------------------------------------------------
!Si se aplica el modelo conjunto y tiene afección la condición similar a la condición seco-mojado, se debe tener en cuenta que las condiciones
!de contorno que se impondrán en el contorno móvil (que no coincidirán con los valores de la solución temporal previa y que son 
!generadas por el modelo subterráneo) se corresponden con un tiempo igual al tiempo para el que se busca la solución.
if (sino.eq.'si') then
call timeasu (i,j,x,y,At)
!Para formulación BG las matrices M (timeasu) valen para las ecuaciones de aguas someras y de N-S 2D. 
!Para formulación estabilizada esto no es así para las nuevas matrices de masa estabilizadas porque se usan parámetros distintos. 
!De todos modos no funciona bien la consideración de matrices de masa estabilizadas y de momento no se usan (líneas comentadas debajo).
 if (navier.eq.'si') then 
  !if (est.eq.'si')then
  !call timeasupgns (i,j,x,y,vv,At)
  !else
  continue
  !endif
 else
  call timeasd (i,j,x,y,At)
  !if (est.eq.'si')then
  !call timeasupgas (i,j,x,y,vv,At)
  !endif
 endif
!Se guarda la matriz de rigidez con las matrices de masa.
do u=1,ndim
cia(u)=ita(u)
ca(u)=sa(u)
enddo
!Se ordenan los coeficientes (alguno podría ser nulo), se ensamblan y se eliminan los coeficientes nulos. Resumiendo, se da el formato MSR.
!Sólo así es posible calcular el producto de la matriz de rigidez (con las matrices de masa) por la solución en el instante anterior.
call orden (inc,k)
do u=1,3*i
 do uu=ita(u),ita(u+1)-1
 voc(u)=voc(u)+sa(uu)*vib(ita(uu))
 enddo
 voc(u)=voc(u)+sa(u)*vib(u)
enddo
!Se vuelve a tener la matriz de rigidez en 'ita,sa' sin ordenar. Así, será posible sumar directamente los coeficientes correspondientes a otras
!matrices que estén en la misma caja de la matriz (por ejemplo A y M) sin buscar donde deben ser almacernardos (dado que el orden con que se toman 
!los elementos de la malla será el mismo). 
do u=1,ndim
ita(u)=cia(u)
sa(u)=ca(u)
enddo
endif 	   

!Se calculan las matrices A y B de las ecuaciones dinámicas. Valen para las ecuaciones de aguas someras y de N-S 2D:
!-------------------------------------------------------------------------------------------------------------------
call cajasab (i,j,x,y,nu,del)  

!En la primera iteración se sustituyen las CC en el vector inicial para calcular las integrales de contorno.
!-----------------------------------------------------------------------------------------------------------
!En el resto de iteraciones se utilizará la solución y los valores serán los mismos a los de la CC al haber reducido el tamaño del sistema con ellas.
if ((it.eq.0).and.(itt.eq.0).and.(navier.eq.'si'))then
!En caso de (navier.eq.'no') esto ya se ha hecho en la subrutina mallaaguassomeras. Se recuerda que en este caso se hace la selección inicial 
!y es necesaria esa subrutina. Además, en este caso no sería necesario hacer esto si la CI o aproximación inicial tiene los valores de las CC.
 do u=1,3*i
  if (vo(u).ne.sqrt(2.0_8)) then 
   if (u.gt.2*i)then
   vv(u)=vo(u)-z(u-2*i)
   vb(u)=vo(u)
   else
   vv(u)=vo(u)
   endif  
  endif
 enddo
endif
 
!Otras matrices y vectores para las ecuaciones de Navier-Stokes 2D:
!------------------------------------------------------------------
if (navier.eq.'si') then        
 !Cálculo de las integrales de contorno	para el caso implícito.
 call vectorcontornopresiones(i,j,x,y,vv,vb,veac,veoc,veuc,nu,ten) 
 if (imp.eq.'no')then
 !Suma del cálculo de otras integrales de contorno para el caso semi-implícito.
 call vectorcontornopresiones(i,j,x,y,vit,vib,veac,veoc,veuc,nu,ten)
 endif
  !Cálculo de las matrices Bt de la ecuación de continuidad
  call cajasbt (i,j,x,y,del)
  !Suma de las integrales de contorno al término independiente. 
  do u=1,i
  vec(u)=del*voc(u)+veac(u)/del
  vec(i+u)=del*voc(i+u)+veoc(u)/del
  vec(2*i+u)=voc(2*i+u)
  enddo
   !Suma de una matriz (más completa) por un vector al término independiente para el caso semi-implícito.
   if (imp.eq.'no')then
    do u=1,ndim
    cia(u)=ita(u)
    ca(u)=sa(u)
    enddo
    call matriznolineal(i,j,x,y,vit,del)
    if (est.eq.'si')then
    call cajasupgns(i,j,x,y,vit,nu,del)
    endif
    call orden (inc,k)
    do u=1,3*i
     do uu=ita(u),ita(u+1)-1
     vec(u)=vec(u)-sa(uu)*vib(ita(uu))
     enddo
     vec(u)=vec(u)-sa(u)*vib(u)
    enddo
    do u=1,ndim
    ita(u)=cia(u)
    sa(u)=ca(u)
    enddo
   endif
  !Cálculo de las matrices estabilizadas.
  if (est.eq.'si')then
  call cajasupgns(i,j,x,y,vv,nu,del)
  endif    

!Otras matrices y vectores para las ecuaciones de aguas someras:
!---------------------------------------------------------------
else
 !Cálculo de las integrales de contorno	y las integrales de la pendiente de fricción para el caso implícito.
 call vectorcontornopresiones(i,j,x,y,vv,vb,veac,veoc,veuc,nu,ten)
 call f(i,j,x,y,z,ma,yn,vv,fx,fy)
  !Incluyen las integrales de la pendiente de fricción estabilizadas.
  if (est.eq.'si')then
  call fsupg(i,j,x,y,ma,vv,fx,fy,fz)
  endif
 if (imp.eq.'no')then
 !Suma del cálculo de otras integrales de contorno y otras integrales de la pendiente de fricción para el caso semi-implícito.
 call vectorcontornopresiones(i,j,x,y,vit,vib,veac,veoc,veuc,nu,ten)
 call f(i,j,x,y,z,ma,yn,vit,fx,fy)
  if (est.eq.'si')then
  call fsupg(i,j,x,y,ma,vit,fx,fy,fz)
  endif
 endif
 !Cálculo de las integrales correspondientes a la lluvia.
 call vectorcontornofuente(i,j,x,y,vef,ql)	
  !Suma de las integrales ánteriores al término independiente.
  do u=1,i
  vec(u)=del*voc(u)-fx(u)*9.81_8/del+veac(u)/del
  vec(i+u)=del*voc(i+u)-fy(u)*9.81_8/del+veoc(u)/del
  vec(2*i+u)=del*voc(2*i+u)-fz(u)*9.81_8/del-veuc(u)/del !+vef(u) !es
  enddo  
   !Suma de una matriz (más completa) por un vector al término independiente para el caso semi-implícito.
   if (imp.eq.'no')then
    do u=1,ndim
    cia(u)=ita(u)
    ca(u)=sa(u)
    enddo
    call cajasde (i,j,x,y,vit,del)
    call matriznolineal(i,j,x,y,vit,del)
    if (est.eq.'si')then
    call cajasupgas(i,j,x,y,vit,nu,del)
    endif
    call orden (inc,k)
    do u=1,3*i
     do uu=ita(u),ita(u+1)-1
     vec(u)=vec(u)-sa(uu)*vib(ita(uu))
     enddo
     vec(u)=vec(u)-sa(u)*vib(u)
    enddo
    do u=1,ndim
    ita(u)=cia(u)
    sa(u)=ca(u)
    enddo
   endif  
  !Cálculo de las matrices D y E de la ecuación de continuidad.
  call cajasde (i,j,x,y,vv,del)
  !Cálculo de las matrices estabilizadas.
  if (est.eq.'si')then
  call cajasupgas(i,j,x,y,vv,nu,del)
  endif  
endif

!Se calculan las matrices C de las ecuaciones dinámicas. Valen para las ecuaciones de aguas someras y de N-S 2D.
!---------------------------------------------------------------------------------------------------------------
call matriznolineal(i,j,x,y,vv,del)	  

!Se ha aplicado un esquema semi-implícito determinado (casi de Crank-Nicolson)
!Con vit se tiene otra versión según el libro de Quarteroni y Valli (chap. 13). De programarla, sólo llamaría una vez a cada caja con vit 
!y me traería hasta aquí el trozo "donde se suma una matriz (mas completa) por un vector al término independiente". 
!En caso de hacerlo así, no podría tener programado a la vez el implícito.

!Método de Picard:
!-----------------
if (((it.lt.5).and.(newton.eq.'si')).or.(newton.eq.'no')) then
 !Se imponen las CC sobre el sistema y se forma la matriz de dimensión (2*i + nºnodos esquina)-nº de CC para en el dominio;
 !Se reduce también en los nodos donde eval=sqrt(3)	si la condición seco-mojado (o la condición similar a la condición seco-mojado) tiene afección.
 call reducciondelsistema(i,c,vec,vu,vn,vb,inc,bcg,ndim)	  												
 !Resolucion mediante gradientes biconjugados. Se utilizará la solución anterior (nonzero=.true.) para acelerar el método.																		
 allocate (vbdin(c),vecdin(c))	 
 do u=1,c			
 vecdin(u)=vec(u)
  if (nonzero)then
  vbdin(u)=vb(u)  
  else
  vbdin(u)=0.0_8
  endif
 enddo
  !Uso de las subrutinas (señaladas en Press et al.) que calculan con un precondicionador diagonal.
  if (bcg.eq.'no') then
  call gradientesbiconjugados(vecdin,c,vbdin,nonzero) 
  !Uso de las subrutinas (incluidas en la librería SLATEC) que calculan con un precondicionador LU aproximado.
  !Allí se definen 'nu' como el nº de coeficientes no nulos en la parte triangular superior de la matriz de rigidez, 'nl' el nº en la parte inferior 
  !y 'nelt' el número total de coeficientes no nulos. Se usan las dimensiones nl+nu+8*n para el vector de coeficientes enteros y nl+nu+4*n+2 para el 
  !vector de coeficientes reales. Se cumplirá que nu+nl = nelt-n = ndim-1-n = ndim-1-c. 
  else                                  	
  call dslubc (vecdin,c,vbdin,ndim+7*c,ndim+3*c+12)	  
  endif
  do u=1,c			
  vb(u)=vbdin(u)
  enddo
 !Arreglo del vector solución introduciendo los valores prescritos (que incluye valores nulos de velocidad y alturas para la parte de la malla seca 
 !en nodos cuadráticos y lineales si tiene afección la condición seco-mojado o la condición similar a la condición seco-mojado) 
 !y valores nulos de altura en los nodos cuadráticos (sólo en los del trozo mojado de la malla en caso de afección de la condición).	
 do u=1,3*i
  if (vu(u).ne.sqrt(2.0_8)) then 
  !Si la condición (vu.ne.) con mayor u es para u<<3*i se entra con c=3*i-1. Se sale con c=3*i y al bucle 'u' aún le quedarán iteraciones (no se entrará).
  !si la condición (vu.ne.) con mayor u es para u=3*i-1 se entra con c=3*i-1. Se sale con c=3*i y queda una última iteración de 'u' (no se entrará).
  !Si la condición (vu.ne.) con mayor u es para u=3*i	se entra con c=3*i-1. Así, no se opera el siguiente if, se saldrá con c=3*i.
  !Al ir por orden 'u' sólo puede ser una unidad mayor a 'c' y 'vb' ya estará ordenado.    
   if (u.le.c) then
   do uu=c,u,-1
   vb(uu+1)=vb(uu)			  
   enddo
   endif
   c=c+1
   vb(u)=vu(u)
  endif 
 enddo

!Método de Newton: sólo se puede usar para formulación BG, para ec. aguas someras, para esquema implícito y para it>4
!--------------------------------------------------------------------------------------------------------------------
else
 !Aplicación del método a la ec. de aguas someras del mismo modo que se aplica a la ec. de N-S 2D en Reddy y Gartling. 
 !Suma de la matriz de rigidez por la solución anterior del problema (rb=vb) y del antiguo término independiente cambiado de signo 
 !al nuevo término independiente.
 do u=1,ndim
 cia(u)=ita(u)
 ca(u)=sa(u)
 enddo
 call orden (inc,k)
 do u=1,3*i
  do uu=ita(u),ita(u+1)-1
  ru(u)=ru(u)+sa(uu)*rb(ita(uu))	   
  enddo
  ru(u)=ru(u)+sa(u)*rb(u)-vec(u)
 enddo
 do u=1,ndim	
 ita(u)=cia(u)
 sa(u)=ca(u)
 enddo 
 !Se calcula la matriz de las derivadas (matriz jacobiana), la nueva matriz del sistema.
 !Dentro de la subrutina se siguen los mismos pasos que en la subrutina reducciondelsistema, sólo que en este caso se 
 !imponen CC nulas allí donde se conozca una CC.
 call jacob (i,j,x,y,vv,nu,ma,ru,vn,vu,vb,c,ten,bcg,ndim)
 allocate (vbdin(c),vecdin(c))	 
 !Resolucion mediante gradientes biconjugados. 
 !Mejor comportamiento inicializando a cero (nonzero=.false.) que guardar la última solución (la solución se hara cero en la convergencia). 
 !Además se evita almacenar una nueva variable ya que en vb no está la última solución del sistema, (ha sido 
 !sobreescrito en la iteración previa con vb=rb-vb y tiene la solución anterior del problema).
 nonzero=.false.
 do u=1,c			
 vecdin(u)=ru(u)  
 vbdin(u)=0.0_8 
 enddo
  if (bcg.eq.'no') then
  call gradientesbiconjugados(vecdin,c,vbdin,nonzero)
  else 
  call dslubc (vecdin,c,vbdin,ndim+7*c,ndim+3*c+12)	   
  endif
  !La solución del sistema tendrá todos los valores más próximos a cero tanto mayor número de iteraciones se hagan.
  do u=1,c			
  vb(u)=vbdin(u)
  enddo
 !Arreglo del vector solución del mismo modo que al terminar cada iteración de Picard (los valores prescritos ahora también serán valores nulos).
 do u=1,3*i
  if (vu(u).ne.sqrt(2.0_8)) then
   if (u.le.c) then
   do uu=c,u,-1
   vb(uu+1)=vb(uu)
   enddo
   endif
   c=c+1
   vb(u)=0.0_8
  endif 
 enddo
 !Nueva solución = solución anterior - solución del sistema
 do u=1,3*i
 vb(u)=rb(u)-vb(u)
 enddo
endif

!Postproceso tras la resolución del sistema:
!-------------------------------------------
!Se calcula el calado en los nodos cuadráticos.
!Da igual que los calados de la solución de aguas someras sean menor que cero porque 
!con la condición seco-mojado se creará una nueva malla que no considerará estos valores.
do u=1,3*i
vt(u)=vb(u)
enddo
do u=1,i	  
if ((vnin(u).ne.1).and.(eval(u).eq.0.0)) then
vt(2*i+u)=vt(2*i+u)-z(u)		                         		   
endif
enddo
it=it+1

!Se muestra la solución de velocidades y calados por pantalla (sólo en los nodos del dominio superficial si la condición
!seco-mojado o la condición similar a la condición seco-mojado tiene afección).	Impresión de resultados fraccionada, de 200 en 200
!líneas (necesario definir 'ua' como entero).
!write(6,*)' '
!ua=1
!uu=0
!do u=1,i
! if (eval(u).eq.0.0) then								 
! write(6,30)'u',u,vt(u),'v',u,vt(u+i)	
! uu=uu+1
!  if (uu.eq.200*ua) then
!  write(6,*)'Mostrados arriba 200 resultados. Pulsa enter.'
!  read(5,*)
!  ua=ua+1
!  endif  
! endif
!enddo
!write(6,*)' '
!ua=1
!uu=0
!do u=1,i
! if ((vnin(u).ne.1).and.(eval(u).eq.0.0)) then		             
! write(6,30)'h',u,vt(u+2*i)
! uu=uu+1
!  if (uu.eq.200*ua) then
!  write(6,*)'Mostrados arriba 200 resultados. Pulsa enter.'
!  read(5,*)
!  ua=ua+1
!  endif
! endif
!enddo
	 
!Se muestra por pantalla cuándo se da el caso particular de tener la solución del flujo de Stokes 2D, en el que no hay convección (se da ya que
!la matriz convectiva no lineal es nula al introducir velocidades nulas y no es necesario resolver un problema no lineal). 
!if ((it.eq.1).and.(itt.eq.0).and.(navier.eq.'si')) then 
!write(6,*)' '
!write(6,*)'Primera iteracion del metodo de Picard con ecuaciones de NS 2D partiendo'
!write(6,*)'de un valor nulo. Equivale a la solucion del flujo de Stokes 2D.'	  
!endif
	  		  
!Control para la parada del bucle para la iteración no lineal.
write(6,*)'Error entre iteraciones no lineales: ',maxval(abs(vv-vt))
!Cuando la diferencia entre cualquier componente del vector solución y la obtenida en la iteración anterior sea menor que la variable tol
!definida en la subrutina aguassomerassubt se detendrán las iteraciones (se sobreescribe aquí tol y así no se entrará en esta subrutina).
if (maxval(abs(vv-vt)).lt.tol) then
tol=1.0
endif

if (bcg.eq.'no') then
deallocate(vbdin,vecdin,ita,sa)
else
deallocate(vbdin,vecdin,cia,ca,cja) 
endif

!Condición seco-mojado ó condición similar a la condición seco-mojado:
!---------------------------------------------------------------------
!La condición seco-mojado deja la solución preparada para otra iteración del modelo superficial (da valores, selecciona dominio y da CC en contorno móvil)
!Es para el modelo superficial.
!La condición similar deja la solución preparada para el modelo subterráneo (da valores). Antes de resolver el subterráneo se usa una subrutina previa
!(que selecciona dominio y da CC). Tras resolver el subterráneo también será necesaria otra subrutina previa (que seleccione dominio y de CC). 
!Es para el modelo conjunto.
do u=1,i
vth(u)=vt(2*i+u)
enddo
!Si se busca la solución de las ecuaciones de N-S 2D no se aplica esta condición. Se aplica si se busca la solución de las ec. de aguas someras
!con o sin iteraciones previas de N-S 2D (siempre sucede al aplicar el modelo conjunto). 
!Por tanto siempre que se entra aquí, se habrá evaluado la condición inicial o la aproximación inicial en la subrutina aguassomerassubt 
if ((navier.eq.'no').or.((navier.eq.'si').and.(ap.eq.'si'))) then 
if (maxval(vth).le.0.0) then 
 !Para it=1 (ya es 1 en la tras la 1ª iteración) si se resuelve con toda la malla y salen todos los calados negativos se tendrá maxval(vth).lt.0.0.
 !Para it>1 se tendrá un dominio menor (con vth=0 fuera), y si en él salen todos los calados negativos se tendrá maxval(vth).le.0.0.   
 write(6,*)'El flujo superficial ha desaparecido'
 do u=1,i
 eval(u)=sqrt(3.0_8)
 enddo
 if (modelo.eq.'superficial') then
  do u=1,3*i
  vt(u)=0.0_8
  enddo
  tol=1.0
 endif
elseif ((minval(vth).lt.0.0).or.(maxval(eval).ne.0.0))then
 !Si se ha considerado toda la malla y los valores de calado son positivos no se entra aquí. 
 !El utilizar una condición que evalúe el dominio considerado (eval) además de otra que evalúe el signo del calado (vth), permite entrar y añadir 
 !elementos al nuevo dominio aún cuando no se obtengan calados negativos (cuando el dominio superficial es menor al dominio de toda la malla).   
 !'eval' se ha evaluado inicialmente en la subrutinas aguassomerassubt
 !>Se aplica la condición: 
 !>Modelo superficial - velocidad nula en contorno móvil 
 !'eval' sólo se ha sobreescrito en la subrutina mallaaguassomeras para it=1.
 !Si algún coeficiente no es nulo se entra aquí.
 !Si todos sus coeficientes son nulos (se ha considerado todo el dominio) y en alguna iteración hay valores negativos de calado negativo,
 !se entrará aquí y 'eval' se modificará dentro de la subrutina. 
 !Una vez que se entra siempre se entrará a no ser que en nuevamalla se den valores nulos a todos los coeficientes.
 !>Modelo conjunto - no se aplica velocidad (posteriormente, en la subrutina mallaaguassomeras, se dan velocidades subterráneas) en contorno móvil 
 !Para it=1 se habrán aplicado condiciones de velocidad nula (depende de la condición inicial o aproximación inicial) en el contorno móvil
 !En it=2 ó it>2 del modelo conjunto siempre se habrán aplicado valores de velocidad subterránea calculados con el modelo subterráneo.
 !Aunque la malla sea mayor o menor se habrán aplicado en todos los nodos del contorno móvil al hacer 1 it por cada convergencia del modelo subterráneo. 
 !'eval' se sobreescribe en la subrutina mallaaguassomeras para cada iteración.
 !Si algún coeficiente no es nulo se entra aquí.
 !Si todos sus coeficientes son nulos (se ha considerado todo el dominio) y hay valores negativos de calado negativo,
 !eval se modificará en la subrutina mallaaguassomeras.		  
 call nuevamalla(modelo,i,j,z,vt,vth,vb,eval,vo)
 !Así, 'vo' seguramente será diferente la proxima vez que se resuelva una iteración.		
endif
endif
deallocate(vn,vv,vec,voc,vu,vef,fx,fy,fz,vth,veac,veoc,veuc,ru,rb)	  
end

!----------------------------------------------------------------------------------------------------------------------------------------------------------
!Subrutina AGUASSUBTERRANEA (1 ecuación de continuidad con las características de conductividad, ángulo de anisotropía y porosidad). 
!En esta subrutina se calcula el nivel freático (altura de la lámina de agua). 
!Se contruye el sistema con la matriz de rigidez (una caja). Se calculan los niveles freáticos en los nodos esquina. 
!Las velocidades de Darcy se calculan postproceso. 
!Como solución en cada iteración habrá altura subterránea y una velocidad (en vtv=vbv y velx,vely).
!Se harán todas las iteraciones no lineales para la convergencia del modelo.
!'it' es la iteración del método para resolver la no linealidad del sistema (sólo se usa Picard). 
!'tol' gestiona la parada del bucle cuando dos soluciones consecutivas son similares.
!'vnv,vov,vuv' son contadores de nodos cuadráticos, lineales o con condición de contorno
!'veic' es el vector que contiene las integrales de contorno, 'vefv' es el término fuente de caudal por lluvia para flujo subterráneo
!'qxx,qyy' llevan los valores de condiciones de contorno de caudal en los nodos donde éstas se conocen (guardadas en qx,qy), 
!y valores de caudales por metro lineal calculados a través de las velocidades subterráneas en el resto de nodos.
!'vecv,vocv' forman el término independiente del sistema  
!'vvv,vtv,vbv' son los vectores donde se guarda la solución del sistema en cada iteración, 'c' es el tamaño del sistema en cada iteración.
!'inc' es la dimensión de la matriz del sistema cuadrada si se usasen elementos cuadráticos para niveles freáticos (nºecuaciones*i)
!Otras variables con el mismo nombre son comentadas en la subrutina aguassomeras.                                                                              
!---------------------------------------------------------------------------------------------------------------------------------------------------------
subroutine aguassubterranea(i,j,x,y,zp,nd,kix,kiy,ag,sino,bcg,At,vnin,vov,eval,ql,qb,qx,qy,vibv,vtv,velx,vely) 
use allocatacion
integer*4, dimension(:),allocatable::vnv
integer*4 i,j,u,uu,c,it,vnin(i),inc,ndim,k			 
real*8, dimension(:),allocatable::vbdin,vecdin,vvv,vefv,vuv,vecv,vocv,vbv,qxx,qyy,veic		  	 
real*8 x(i),y(i),zp(i),tol,eval(i),At,vov(i),ag(i)  
real*8 vtv(i),vibv(i),kix(i),kiy(i),nd(i),qx(i),qy(i),qb(i),ql(i),velx(i),vely(i)									  		
logical nonzero	  
character sino*2,bcg*2

allocate(vnv(i),vvv(i),vefv(i),vuv(i),vecv(i),vocv(i),vbv(i),qxx(i),qyy(i),veic(i)) 

30  format(2(3x,A2,'(',I5,')=',E15.8E2))	  
	 
!Inizialización previa de variables:
!-----------------------------------
inc=i
do u=1,i
vbv(u)=0.0_8
enddo

!Inicio de iteraciones para un incremento de tiempo o para el problema estacionario (se aplica un esquema implícito):
!-------------------------------------------------------------------------------------------------------------------- 
it=0
tol=1e-6
do while (tol.eq.1e-6) 

!Selección de en qué iteraciones se utilizará la solución anterior para el método de los gradientes biconjungados
!----------------------------------------------------------------------------------------------------------------
if (it.eq.0)then
nonzero=.false.	
else
nonzero=.true.
endif

!Inizialización de variables (procedimientos ya aplicados con la subrutina aguassomeras):
!----------------------------------------------------------------------------------------
 do u=1,i
 vvv(u)=vtv(u)
 vocv(u)=0.0_8
 vefv(u)=0.0_8
 veic(u)=0.0_8
 enddo
 !El vn original es ahora vnin. El nuevo vn será diferente si el dominio es menor al de toda la malla.
 !En la ecuación se eliminarán las ecuaciones relativas a todos los nodos cuadráticos y no sólo los del dominio suterráneo.
 !Con vn se eliminarán los de este dominio y el resto se eliminarán con las CC al no considerar el resto de la malla.
 do u=1,i
  if ((vnin(u).eq.1).and.(eval(u).eq.0.0)) then
  vnv(u)=1
  else
  vnv(u)=0
  endif
 enddo
 !Se dimensionan ita y sa (también cia y ca), de acuerdo al número total de coeficientes que se generarán con las matrices elementales.
 !Además ita ya lleva referenciados el número de coeficientes que habrá por fila (en sus primeras i+1 componentes) antes de dar el formato MSR.
 !Cuando se calcule una matriz se almacenarán sus coeficientes directamente en estos vectores.
 call dimvectsb (i,j,vnin)
 ndim=ita(inc+1)-1
	
!En caso de resolver de forma transitoria con incrementos de tiempo (esquema implícito considerando las matrices de masa).
!-------------------------------------------------------------------------------------------------------------------------
!Si se aplica el modelo conjunto y tiene afección la condición similar a la condición seco-mojado, se debe tener en cuenta que las condiciones
!de contorno que se impondrán en el contorno móvil (que no coincidirán con los valores de la solución temporal previa y que son 
!generadas por el modelo superficial) se corresponden con un tiempo igual al tiempo para el que se busca la solución.
if (sino.eq.'si') then
call timesubt (i,j,x,y,nd,At)
!Otra opción es utilizar la siguiente subrutina en vez de la anterior, calculando la matriz de masa concentrada:
!call timesubtconc(i,j,x,y,nd,At,ndim) !cñ
do u=1,ndim
cia(u)=ita(u)
ca(u)=sa(u)
enddo
call orden (inc,k)
do u=1,i
 do uu=ita(u),ita(u+1)-1
 vocv(u)=vocv(u)+sa(uu)*vibv(ita(uu))
 enddo
 vocv(u)=vocv(u)+sa(u)*vibv(u)
enddo
do u=1,ndim
ita(u)=cia(u)
sa(u)=ca(u)
enddo
endif
do u=1,i
vtv(u)=0.0_8
enddo						

!Arreglo de las condiciones de contorno.
!---------------------------------------
!No es necesario sustituir las CC en el vector inicial ya que las integrales de contorno no dependen de las variables del sistema.
!Para el cálculo de las integrales de contorno se usan los valores de caudal qx,qy (en los nodos lineales) leídos inicialmente del fichero. 
!En vez de dar valores nulos donde no hay condición (por ejemplo poniendo qxx=0, qyy=0 si qx(u),qy(u)=sqrt(3)), se calculan unos valores aproximados de 
!caudal (con valores de velocidades calculados post-proceso en la iteración anterior). No se tendrían a través de la solución ya que esta variable no es 
!solución del sistema. Sin embargo, éstos no tendrán influencia alguna, pues los coeficientes calculados con ellos desaparecen por ensamblar 
!vectores elementales o por aplicar condiciones de nivel freático sobre el sistema (en el contorno móvil siempre se aplican). 	 
do u=1,i		   
 if (qx(u).eq.sqrt(3.0_8))then	   
  !Se tiene en cuenta la condición de espesor mínimo que también aplica en la subrutina cajasasubt.
  !De este modo se evita utilizar espesores negativos cuando aparecen (teniendo esta condición un fin parecido al que se persigue con la aplicación
  !de la condición seco-mojado para el modelo superficial).				
  if ((vvv(u)-zp(u)).le.0.0) then
  qxx(u)=velx(u)*0.001_8
  else
  qxx(u)=velx(u)*(vvv(u)-zp(u))
  endif
 else
 qxx(u)=qx(u)
 endif
 if (qy(u).eq.sqrt(3.0_8))then
  if ((vvv(u)-zp(u)).le.0.0) then
  qyy(u)=vely(u)*0.001_8
  else
  qyy(u)=vely(u)*(vvv(u)-zp(u))
  endif
 else
 qyy(u)=qy(u)
 endif
enddo

!Cálculo de otras matrices y vectores para la ecuaciones subterránea:
!--------------------------------------------------------------------
!En la siguiente subrutina se calculan las integrales de contorno.
call vectorcontornocaudales(i,j,x,y,veic,qxx,qyy)
!Cálculo de las integrales correspondientes a la lluvia.
call vectorcontornofuentesub(i,j,x,y,vefv,ql)
do u=1,i
vecv(u)=vocv(u)-veic(u)-qb(u)+vefv(u)
enddo
!Cálculo de la matriz A (Asx+Asy).
call cajasasubt (i,j,x,y,zp,vvv,kix,kiy,ag)

!Se inicializa vuv, para imponer las condiciones de contorno sobre el sistema
!----------------------------------------------------------------------------
do u=1,i
vuv(u)=vov(u)
enddo

!Aplicación del método de Picard (procedimientos ya aplicados con la subrutina aguassomeras):
!--------------------------------------------------------------------------------------------
call reducciondelsistema(i,c,vecv,vuv,vnv,vbv,inc,bcg,ndim)				
!Resolucion mediante gradientes biconjugados
allocate(vbdin(c),vecdin(c))  
do u=1,c			
vecdin(u)=vecv(u)
!Si (nonzero), que se da para it>0, vbv tendrá coeficientes no nulos. En otro caso siempre será nulo. 
vbdin(u)=vbv(u)
enddo
if (bcg.eq.'no') then
call gradientesbiconjugados(vecdin,c,vbdin,nonzero)
else
call dslubc (vecdin,c,vbdin,ndim+7*c,ndim+3*c+12)	  
endif
do u=1,c			
vtv(u)=vbdin(u) 
enddo

!Arreglo del vector solución introduciendo los valores prescritos de hd (y hd nulo en zonas donde hay solución aguas someras
!tanto en nodos cuadráticos como lineales) y los valores nulos en los nodos cuadráticos (de la nueva malla).
!En caso de aplicar el modelo conjunto se tendrán valores nulos de las variables donde no se calcule flujo subterráneo.... si tiene afección la cond
do u=1,i
 if (vuv(u).ne.sqrt(2.0_8)) then
  if (u.le.c) then
  do uu=c,u,-1
  vtv(uu+1)=vtv(uu)		   
  enddo
  endif
  c=c+1
  vtv(u)=vuv(u)
 endif 
enddo

!Postproceso tras la resolución del sistema:
!-------------------------------------------
!Escritura de la solución sobre el vector vbv.
do u=1,i
vbv(u)=vtv(u)
enddo
!Cálculo de las velocidades tras cada iteración.
call velocidadessubterraneas (i,j,x,y,kix,kiy,ag,vtv,velx,vely)
it=it+1

!Se muestra la solución de velocidades subterráneo y niveles freáticos por pantalla (sólo en los nodos del dominio subterráneo si la 
!condición similar a la condición seco-mojado tiene afección).	Impresión de resultados fraccionada, de 200 en 200
!líneas (necesario definir 'ua' como entero).
!write(6,*)' '
!ua=1
!uu=0
!do u=1,i
! if ((vnin(u).ne.1).and.(eval(u).eq.0.0)) then		                     
! write(6,30)'us',u,velx(u),'vs',u,vely(u)
! uu=uu+1
!  if (uu.eq.200*ua) then
!  write(6,*)'Mostrados arriba 200 resultados. Pulsa enter.'
!  read(5,*)
!  ua=ua+1
!  endif
! endif	
!enddo
!write(6,*)' '
!ua=1
!uu=0
!do u=1,i
! if ((vnin(u).ne.1).and.(eval(u).eq.0.0)) then							 
! write(6,30)'hd',u,vtv(u)
! uu=uu+1
!  if (uu.eq.200*ua) then
!  write(6,*)'Mostrados arriba 200 resultados. Pulsa enter.'
!  read(5,*)
!  ua=ua+1
!  endif
! endif
!enddo

!Control para la parada del bucle para la iteración no lineal.
write(6,*)'Error entre iteraciones: ',maxval(abs(vvv-vtv))
!Cuando la diferencia entre cualquier componente del vector solución y la obtenida en la iteración anterior sea menor que la variable tol
!se detendrán las iteraciones (se sobreescribe aquí tol y así se saldrá del bucle dentro de esta subrutina).
if (maxval(abs(vvv-vtv)).lt.tol) then
tol=1.0
write(6,*)'Convergencia modelo subterraneo en iteracion:', it
endif

if (bcg.eq.'no') then
deallocate(vbdin,vecdin,ita,sa)
else
deallocate(vbdin,vecdin,cia,ca,cja)  
endif 
enddo

!Si se aplica el modelo conjunto y tuvo afección la condición similar a la seco-mojado se habrán dado condiciones de contorno
!en el contorno móvil (previamente, a través de la subrutina mallasubterranea, se dan niveles freáticos).
!Así, 'vov' seguramente será diferente la proxima vez que se resuelva esta ecuación. 
								  
deallocate(vnv,vvv,vefv,vuv,vecv,vocv,vbv,qxx,qyy,veic)
end

!----------------------------------------------------------------------------------------------------------------------------------------------
!Subrutinas de gestión del movimiento de contornos móviles.
!----------------------------------------------------------------------------------------------------------------------------------------------
!Subrutina MALLAAGUASSOMERAS
!Esta subrutina puede modificar la malla del dominio superficial (sobreescribiendo el fichero malla.txt) haciendo una selección del dominio 
!superficial en base al dominio subterráneo, localizando el mismo contorno móvil definido con la subrutina mallasubterranea.
!Además da CC de velocidad subterránea (o nula si no la hay en la CI o aproximación inicial o no se ha calculado solución subterránea) 
!en el contorno móvil.
!Se aplica antes de resolver/aplicar la ec superficial.
!Se utiliza una vez (1 it) si se aplica el modelo superficial para definir inicialmente el dominio superficial 
!dando CC de valor nulo. Las siguientes selecciones (por iteración) del dominio superficial corresponden a la subrutina nuevamalla.
!Se utiliza en cada iteración (en cada iteración conjunta) si se aplica el modelo cunjunto (en la primera iteración también se define el 
!dominio y se dan CC de valor nulo). 
!----------------------------------------------------------------------------------------------------------------------------------------------
subroutine mallaaguassomeras(i,j,z,vta,velx,vely,eval,vuc,vt)
use interaccion
integer*4 i,j,u,ui,uuu,uu,io		   
real*8 z(i),vt(3*i),vta(i),velx(i),vely(i),eval(i),vuc(3*i),fi												   
character fe*1

10  format(I5)
21  format(3/,I5)								
23  format(4/,A80)									
26  format(6X,6(X,I5))	
27  format(I5,X,6(X,I5))	
46  format(X,I5,X,A1,X,F11.7)						

!Inizialización previa de variables:
!-----------------------------------
io=0
uu=0							  
do u=1,i
eval(u)=sqrt(3.0_8)
enddo
do u=1,3*i		  
vuc(u)=sqrt(2.0_8)	    				   
enddo

!Selección de dominio considerando los valores de nivel freático:
!----------------------------------------------------------------
open(unit=4,file='C:\mallainicial.txt',status='old')
!La opción scratch crea un archivo temporal que desaparece al cerrar el archivo
open(2,status= 'scratch') 
read(4,21)j
read(4,'(A)')a
 do u=1,j
 read(4,26) no(1),no(4),no(2),no(5),no(3),no(6)
   !Selección inicial. Aquí sólo se consideran los elementos con todos sus nodos de nivel freático nulo.   
   if ((vta(no(1)).eq.0.0).and.(vta(no(2)).eq.0.0).and.(vta(no(3)).eq.0.0)) then
   uu=uu+1
   write(2,27)uu,no(1),no(4),no(2),no(5),no(3),no(6)
   !Se considerará el elemento.
   do uuu=1,6
   eval(no(uuu))=0.0_8
   enddo
   endif
 enddo
 rewind(4)
 read(4,23)a
 do u=1,j
 read(4,26)no(1),no(4),no(2),no(5),no(3),no(6)
   do ui=1,3
   b=ui+1-sb(ui)
   c=ui+2-sb(ui+1)
   !Se toman los elementos que rodean a la selección inicial. También se toman elementos pegados a nodos aislados o a líneas que no forman elementos 
   !(regueros) donde no se calculó solución subterránea (niveles freáticos nulos). Es probable encontrar regueros pegados a la selección inicial. 
   !Estos casos no se darán debido a una previa búsqueda de flujo superficial aislado en la solución subterránea, ya que esto no se hace. 
   !Dentro de los elementos que se toman, no se ha calculado el flujo subterráneo (hay solución para el flujo superficial, de la anterior 
   !iteración de la ecuación superficial).   	   
   if ((vta(no(ui)).ne.0.0_8).and.(vta(no(b)).eq.0.0_8).and.(vta(no(c)).eq.0.0_8)) then  
	!Elementos con un nodo apoyado en el contorno móvil - tienen nodos esquina con: dos niveles freático nulos y uno no nulo.
	uu=uu+1
    write(2,27)uu,no(1),no(4),no(2),no(5),no(3),no(6)
	do uuu=1,6
    eval(no(uuu))=0.0_8
    enddo		 	
   elseif ((vta(no(ui)).ne.0.0_8).and.(vta(no(b)).ne.0.0_8).and.(vta(no(c)).eq.0.0)) then
	!Elementos con dos nodos apoyados en el contorno móvil - tienen nodos esquina con: un nivel freático nulo y dos no nulos. 
	!A través de ellos se dan las CC de velocidad subterránea. No hay que considerar elementos con dos niveles freáticos nulos y uno no nulo 
	!para ello pues aún en caso de que haya dos o más elementos juntos de este tipo compartirán este nodo de nivel freático no nulo, y este nodo 
	!formará finalmente parte de un elemento de los analizados aquí.
	uu=uu+1
    write(2,27)uu,no(1),no(4),no(2),no(5),no(3),no(6)
	do uuu=1,6
    eval(no(uuu))=0.0_8
    enddo
	 !Se llega a una interfaz sumando elementos (considerando también el caso particular, se hace a continuación) a la selección inicial y será la 
	 !interfaz que se calculó en la subrutina mallasubterránea. Por tanto, la interfaz es la misma para flujo de aguassomeras y flujo subterráneo. 
	 !Está situada en nodos donde hay calado, en la orilla del dominio superficial. Por tanto se pasan valores que ya tienen esos nodos.
	 !Serán CC de velocidad subterránea. Se da un valor interpolado en los nodos cuadráticos (en ellos no se ha calculado valor).
	 vuc(no(ui))=velx(no(ui))   		  	 
     vuc(no(b))=velx(no(b))
	 vuc(no(ui+3))=0.5_8*(velx(no(ui))+velx(no(b)))
	 vuc(i+no(ui))=vely(no(ui))   		  	 
     vuc(i+no(b))=vely(no(b))
	 vuc(i+no(ui+3))=0.5_8*(vely(no(ui))+vely(no(b)))
	 !Dado que la interfaz (contorno móvil) es la misma, es posible garantizar la conservación. Las velocidades subterráneas estarán mejor 
	 !aproximadas tanto menor sean los elementos pegados a la interfaz, dado que se interpolan linealmente.
	 !En otro caso (dos contornos móviles) el número de nodos en cada interfaz puede ser diferente, y aunque sea el mismo será complicado 
	 !calcular las componentes en nodos diferentes de forma que el caudal se conserve.	   
   endif
   enddo    
 enddo
 
!Ojo, es posible que el flujo subterráneo haga casi desaparecer al flujo superficial, si la dimensión de éste es pequeña (por ejemplo 
!un río de poco de ancho). En este caso, podría seleccionarse el dominio a partir de un reguero, y este reguero estar dividido en dos 
!partes en alguna zona. En este caso habrá que refinar más la malla.
!Aún en caso de que los elementos considerados pegados a ambos trozos de reguero se toquen formando un dominio contínuo habría problemas
!ya que se dará condición de velocidad nula en el punto donde ambas partes se unen.

!Caso particular:
!----------------
!Es posible tener elementos que no se consideran con el procedimiento anterior y que son necesarios para tomar todos los elementos hasta
!el contorno móvil.
!Estos elementos no estarán pegados a la selección inicial, y cumplen la particularidad de tener sus tres nodos esquina apoyados en el contorno móvil.
!Podrían tampoco estar pegados a la última selección hecha si forman una banda con el ancho de un elemento.  
rewind(4)
read(4,23)a
 do u=1,j
 read(4,26) no(1),no(4),no(2),no(5),no(3),no(6)   
   !Se selecciona el elemento (pegado o no pegado a la última selección hecha antes de entrar aquí).
   if (vta(no(1)).eq.(vt(2*i+no(1))+z(no(1)))) then
   if ((vta(no(2)).eq.(vt(2*i+no(2))+z(no(2)))).and.(vta(no(3)).eq.(vt(2*i+no(3))+z(no(3))))) then
   uu=uu+1
   write(2,27)uu,no(1),no(4),no(2),no(5),no(3),no(6)		
   do uuu=1,6
   eval(no(uuu))=0.0_8
   enddo 
    !Se corrijen las condiciones de contorno. Se consideran diferentes casos en función de si el elemento
	!está o no pegado a la última selección hecha (selección que puede estar continuamente cambiando).
	do ui=1,3
    b=ui+1-sb(ui)
    c=ui+2-sb(ui+1)
    e=ui+4-sb(ui)
    f=ui+5-sb(ui+1) 
	 !Elemento no pegado a la última selección.
	 if ((vuc(no(ui+3)).eq.sqrt(2.0_8)).and.(vuc(no(e)).eq.sqrt(2.0_8)).and.(vuc(no(f)).eq.sqrt(2.0_8))) then
	 vuc(no(ui))=velx(no(ui))
	 vuc(i+no(ui))=vely(no(ui))
	 vuc(no(b))=velx(no(b))
	 vuc(i+no(b))=vely(no(b))
	 vuc(no(c))=velx(no(c))
	 vuc(i+no(c))=vely(no(c))
	 vuc(no(ui+3))=0.5_8*(velx(no(ui))+velx(no(b)))
     vuc(i+no(ui+3))=0.5_8*(vely(no(ui))+vely(no(b)))
	 vuc(no(e))=0.5_8*(velx(no(b))+velx(no(c)))
     vuc(i+no(e))=0.5_8*(vely(no(b))+vely(no(c)))
	 vuc(no(f))=0.5_8*(velx(no(ui))+velx(no(c)))
     vuc(i+no(f))=0.5_8*(vely(no(ui))+vely(no(c)))
	 exit
	 !Elemento pegado a la última selección hecha (pegado por uno de sus lados).
	 elseif ((vuc(no(ui+3)).ne.sqrt(2.0_8)).and.(vuc(no(e)).eq.sqrt(2.0_8)).and.(vuc(no(f)).eq.sqrt(2.0_8))) then
	 vuc(no(c))=velx(no(c))
	 vuc(i+no(c))=vely(no(c))
	 vuc(no(ui+3))=sqrt(2.0_8)
     vuc(i+no(ui+3))=sqrt(2.0_8)
	 vuc(no(e))=0.5_8*(velx(no(b))+velx(no(c)))
     vuc(i+no(e))=0.5_8*(vely(no(b))+vely(no(c)))
	 vuc(no(f))=0.5_8*(velx(no(ui))+velx(no(c)))
     vuc(i+no(f))=0.5_8*(vely(no(ui))+vely(no(c)))
     exit
	 !Elemento encerrado por la última selección hecha (pegado por dos de sus lados).
	 elseif ((vuc(no(ui+3)).ne.sqrt(2.0_8)).and.(vuc(no(e)).ne.sqrt(2.0_8)).and.(vuc(no(f)).eq.sqrt(2.0_8))) then
	 vuc(no(ui+3))=sqrt(2.0_8)
     vuc(i+no(ui+3))=sqrt(2.0_8)
	 vuc(no(e))=sqrt(2.0_8)
     vuc(i+no(e))=sqrt(2.0_8)
	 vuc(no(f))=0.5_8*(velx(no(ui))+velx(no(c)))
     vuc(i+no(f))=0.5_8*(vely(no(ui))+vely(no(c)))
	 exit
	 endif
    enddo
   endif
   endif
 enddo

!Valores de la iteración anterior:
!---------------------------------
!Se toman como valores de la iteración anterior (variables de la ecuación superficial) los valores de las CC de velocidad dados.
!Se trata valores en nodos del contorno móvil. De esta forma afectarán a las integrales de contorno (y existirá conservación).
!Ojo, esto da problemas en la convergencia (al modificar la variable 'vt') si se usa la primera opción propuesta a continuación.
do u=1,2*i
 if (vuc(u).ne.sqrt(2.0_8)) then		   
 vt(u)=vuc(u)
 endif
enddo

!Arreglo de las condiciones de contorno:
!---------------------------------------
!Fuera de la nueva malla se dan condiciones de contorno nulas.
!El sistema tendrá coeficientes sólo referidos los nodos de la malla que se selecciona, pero el sistema tendrá la dimensión correspondiente 
!a los nodos de toda la malla (no se pierde la referencia del nodo). Con esto también se reducirá el tamaño del sistema para considerar el 
!sistema correspondiente.
do u=1,i
if (eval(u).eq.sqrt(3.0_8)) then
vuc(u)=0.0_8
vuc(i+u)=0.0_8
vuc(2*i+u)=0.0_8
endif
enddo

!Introducción de las CC de flujo superficial que existen en el contorno no móvil (leídas desde fichero). El contorno no móvil será el trozo de 
!contorno que también encierra a la selección de elementos y que es parte del contorno de toda la malla (del dominio conjunto).
do u=1,i+2
 read(4,'(A)')a
enddo
eaf=0
do while (eaf.eq.0)
 read (4,46,iostat=eaf) u,fe,fi
 io=io+1 
!cñ
!Intersección del contorno móvil con el contorno del dominio conjunto (existe CC en el fichero). 
!>Primera opción que permite la utilización exacta de las CC del fichero. En la intersección:
!Se deja el valor de velocidad subterránea únicamente como valor de la iteración anterior (se usan en las integrales de contorno). 
!Se da la condición de contorno existente (de velocidad o de calado, se usan en la reducción del sistema).
!Respecto a la convergencia, en la itersección, si hay CC de calado la velocidad calculada se modifica evitando 
!la convergencia y si hay CC de velocidad la velocidad de la iteración anterior nunca tiene ese valor. 
!El problema se evita modificando la variable vt debajo, fuera de este bucle "do while", y no arriba (si se modifican valores de la iteración anterior, 
!vuc y vt deben coincidir para que haya convergencia). Pero no se utilizarán las velocidades subterráneas en la intersección (ojo con la conservación).
! if (eaf.eq.0) then 
!  if (io.le.i) then
!   !Hay CC de velocidad en fichero, se sobreescribe la condición en caso de que se haya dado condición de velocidad (intersección).
!   !En otro caso vuc también se sobreescribirá (pues vale sqrt(2.0_8) indicando que no hay condición). 
!   if ((vuc(u).ne.0.0_8).and.(fe.ne.'o')) then 
! 	vuc(u)=fi
! 	endif
!  elseif ((io.gt.i).and.(io.le.2*i)) then
!   !Lo mismo pero para la velocidad en dirección y.
!   if ((vuc(i+u).ne.0.0).and.(fe.ne.'o')) then 
! 	vuc(i+u)=fi
! 	endif															   
!  else
! 	!Hay CC de calado en fichero, se elimina la condición en caso de que se haya dado condición de velocidad (intersección).
! 	!Se da el valor sqrt(2.0_8) para indicar que no hay condición en todos los casos, incluyendo al caso anterior (donde se da CC de calado
!   !en el fichero, nunca se da CC de velocidad).
! 	if ((vuc(2*i+u).ne.0.0_8).and.(fe.ne.'o')) then 	 
! 	vuc(u)=sqrt(2.0_8)
! 	vuc(i+u)=sqrt(2.0_8)
! 	vuc(2*i+u)=fi+z(u)
! 	endif
!  endif
! endif 
!>Segunda opción que asegura la convergencia pero sin utilizar exactamente las CC del fichero. En la intersección:
!Si existe CC calado se añade CC de velocidad subterránea (dar sólo la CC velocidad es complicado porque habría que dar CC en nodos cuadráticos 
!no pertenecientes al contorno móvil).
!Si existe CC de velocidad se substituye por la de velocidad subterránea.
 if (eaf.eq.0) then 
  !Si vuc es distinto de cero (o eval es distinto de sqrt(3.0)) el nodo pertenece a la nueva malla.
  !Si vuc es sqrt(2.0_8) el nodo no tiene valor de contorno (nunca lo tendrá si no se está en el contorno). 
  if (io.le.i) then
    !Hay CC de velocidad en el fichero. Se asegura que sólo exista la CC dada de velocidad.
    if ((vuc(u).eq.sqrt(2.0_8)).and.(fe.ne.'o')) then 
 	vuc(u)=fi
 	endif
  elseif ((io.gt.i).and.(io.le.2*i)) then
    !Lo mismo pero para la velocidad en dirección y.
    if ((vuc(i+u).eq.sqrt(2.0_8)).and.(fe.ne.'o')) then 
 	vuc(i+u)=fi
 	endif															   
  else
    !Ocurre en todo el contorno móvil que vuc(2*u+i) es sqrt(2.0_8) al no darse CC de calado a través de la solución subterránea.
 	!Hay CC de calado en el fichero. Se da la CC superficial en la intersección y no se elimina la que se haya dado de velocidad
	if ((vuc(2*i+u).eq.sqrt(2.0_8)).and.(fe.ne.'o')) then 	
 	vuc(2*i+u)=fi+z(u)
 	endif
  endif
 endif
enddo
		
!Escritura de la malla:
!----------------------
!La malla tendrá la i original porque los números de nodo irán hasta i, y la nueva j para construir bien las cajas existentes
!(no se pierde la referencia del nodo).
rewind(2)                                 
close(4)       
open(1,file = 'C:\malla.txt',status='replace')
write(1,'(A)')'nudos'
write(1,10)i
write(1,'(A)')'elementos'
write(1,10)uu
write(1,'(A)')'elemento p1  p2  p3  p4  p5  p6 (antihorarios)'
do u=1,uu
 read(2,'(A)')a
 write(1,'(A)')a
enddo
close(2) 
close(1)
j=uu
end

!----------------------------------------------------------------------------------------------------------------------------------------------
!Subroutina NUEVAMALLA
!Esta subrutina aplica la condición seco-mojado o la condición similar a la condición seco-mojado. 
!Se aplica tras resolver la ec superficial.
!>Se utiliza (condición seco-mojado) en cada iteración si se aplica el modelo superficial para definir el dominio superficial. 
!Se hace la selección o deselección de elementos para el dominio superficial (en malla.txt) en función de los calados calculados con la 
!ecuación correspondiente y se dan nuevos valores de calado como valores de la iteración anterior. Se da CC de valor nulo en el contorno móvil.  
!>Se utiliza en cada iteración si se aplica el modelo conjunto pero sólo aplica valores de calado de la iteración anterior para una posterior 
!selección del dominio subterráneo con la subrutina mallasubterránea. Cualquier selección de elementos no se considera (no va en malla.txt). 
!Posibilidades (en esta iteración, esto es, cada vez que se entra):
!Selección de elementos pegados al contorno móvil actual (u original). 
!Deselección de cualquier elemento dentro del dominio delimitado por el contorno móvil actual.
!Si el calado es positivo únicamente en los nodos de una línea (regero), éstos no se consideran (habría que refinar la malla).
!----------------------------------------------------------------------------------------------------------------------------------------------
subroutine nuevamalla(modelo,i,j,z,vt,vv,vb,eval,vuc)
use interaccion
integer*4, dimension(:),allocatable::pa,pe
integer*4 i,j,u,ui,uuu,uu,io
real*8, dimension(:),allocatable::si								  
real*8 vt(3*i),vv(i),vb(3*i),eval(i),z(i),vuc(3*i),fi,hmin 
character modelo*12,fe*1

allocate(si(3*i),pa(i),pe(i))

10  format(I5)
21  format(3/,I5)
23  format(4/,A80)
26  format(6X,6(X,I5))	
27  format(I5,X,6(X,I5))
46  format(X,I5,X,A1,X,F11.7)
 												   
!Parámetro para evitar inestabilidades:
!--------------------------------------
!Se definen los elementos interfaz como aquellos elementos pegados al contorno móvil y no considerados en el dominio.
!Es posible que la nueva solución conlleve a deseleccionar elementos seleccionados en una iteración anterior. 
!Este proceso puede ser indefinido si la altura del agua coincide aproximadamente con la cota de los nodos de los elementos interfaz. 
!La solución es usar un parametro 'hmin' que asegure que exista una cierta altura para la selección (se limita la adicción de elemetos).   
!El parametro deberá ser mayor tanto menor sea el número de elementos interfaz.	
!hmin=1e-3   
hmin=2.0_8	!es

!Inizialización previa de variables:
!-----------------------------------
io=0
uu=0
do u=1,3*i		  
vuc(u)=sqrt(2.0_8)
si(u)=0.0_8 
enddo
do u=1,i
eval(u)=sqrt(3.0_8)
enddo

!Evalución de la solución obtenida:
!---------------------------------- 			
!Los nuevos valores dados como valores de la iteración anterior siempre serán erróneos (son supuestos) y más en velocidades al ser más variables. 
!Sin embargo es necesario dar valores de la iteración anterior de calado (y altura) en los nuevos nodos dado que al aplicar la ecuación de aguas 
!someras no puede haber valores nulos o negativos. Se dan en nuevo contorno móvil si se añaden elementos (hay valores nulos), y no se dan en él
!si se eliminan elementos (hay valor). 
!No es necesario dar valores de la iteración anterior de velocidad excepto en el nuevo contorno móvil donde tomarán el mismo valor que las CC.
!Por tanto no se darán en el contorno móvil actual ni en nodos cuadráticos entre contorno móvil actual y nuevo si se añaden elementos (hay valores 
!nulos o de velocidad subterránea en el primer caso y valores nulos en el segundo caso), ni se dan en cualquier nodo perteneciente al dominio pero 
!no al nuevo contorno móvil si se eliminan elementos (hay valor). 
open(unit=4,file='C:\mallainicial.txt',status='old') 
open(2,status= 'scratch')
read(4,21)j
read(4,'(A)')a

!Se decide si con el calado calculado se utiliza el siguiente elemento, no se utiliza, o se eliminan elementos de la malla anterior.
!Esto permitirá, por ejemplo, añadir elementos en una parte del contorno, eliminar elementos en otra parte, y no hacer nada en otra parte. 
!Se trabaja con vv que es vt(2*i+u) original ya que vt(2*i+u) se sobreescribirá. 
  do u=1,j
  read(4,26)no(1),no(4),no(2),no(5),no(3),no(6)
  !Primer proceso. Aquí sólo se consideran los elementos con todos sus nodos positivos. O se seleccionan los mismos elementos de la malla anterior 
  !o se seleccionan menos elementos (eliminación de elementos en este caso). 
   if ((vv(no(1)).gt.0.0).and.(vv(no(2)).gt.0.0).and.(vv(no(3)).gt.0.0)) then  
    uu=uu+1
    write(2,27)uu,no(1),no(4),no(2),no(5),no(3),no(6)
	cycle 
   endif 
   !Segundo proceso. Se analizan elementos pegados al contorno móvil original que no pertenecen a la malla anterior. O no se selecciona o se 
   !selecciona (adicción en este caso). 
   do ui=1,3
    b=ui+1-sb(ui)
    c=ui+2-sb(ui+1)
    e=ui+4-sb(ui)
    f=ui+5-sb(ui+1)
	if ((vv(no(ui)).gt.0.0).and.(vv(no(b)).gt.0.0).and.(vv(no(c)).eq.0.0)) then
	!Elementos con dos nodos apoyados en el contorno móvil original. Si la altura en ambos es superior a la cota en el tercero (y 'hmin') se
    !considera el elemento.
    !Se dan valores como valores de la iteración anterior. En el nuevo nodo se suman los valores de altura de los otros dos, a posibles valores ya 
	!tenidos en cuenta en este proceso. Por ello con la variable 'si' se realizará la media de estas alturas en cada nodo posteriormente.
	if (((vv(no(ui))+z(no(ui))).gt.(z(no(c))+hmin)).and.((vv(no(b))+z(no(b))).gt.(z(no(c))+hmin))) then		   	  		 
	uu=uu+1
	write(2,27)uu,no(1),no(4),no(2),no(5),no(3),no(6)
	vb(2*i+no(c))=vb(2*i+no(c))+vv(no(ui))+z(no(ui))+vv(no(b))+z(no(b))	  
	si(2*i+no(c))=si(2*i+no(c))+2.0_8	
	exit
	endif
    !Elementos con un nodo apoyado en el contorno móvil original. Si la altura en él es superior a la cota en los otros dos nodos del elemento (y 'hmin')
	!se considera el elemento.
	!Se dan valores como valores de la iteración anterior. En los dos nuevos nodos se suman el valor de altura del otro, a posibles valores ya tenidos 
	!en cuenta en este proceso. Por ello con la variable 'si' se realizará la media en cada nodo posteriormente.   
    elseif ((vv(no(ui)).gt.0.0).and.(vv(no(b)).eq.0.0).and.(vv(no(c)).eq.0.0)) then
    if (((vv(no(ui))+z(no(ui))).gt.(z(no(b))+hmin)).and.((vv(no(ui))+z(no(ui))).gt.(z(no(c))+hmin))) then		 
	uu=uu+1
	write(2,27)uu,no(1),no(4),no(2),no(5),no(3),no(6)
	vb(2*i+no(b))=vb(2*i+no(b))+vv(no(ui))+z(no(ui))  
	vb(2*i+no(c))=vb(2*i+no(c))+vv(no(ui))+z(no(ui))  	 
	si(2*i+no(b))=si(2*i+no(b))+1.0_8
	si(2*i+no(c))=si(2*i+no(c))+1.0_8		
	exit
	endif						   	
	endif
   enddo
 enddo
 !Cálculo de valores medios:
 do u=1,i
 if (si(2*i+u).ne.0.0_8) then
 vb(2*i+u)=vb(2*i+u)/si(2*i+u)
 vt(2*i+u)=vb(2*i+u)-z(u)
 endif
 enddo
 !A continuación se sobreescribe la variable 'eval' para que indique cuáles son los elementos considerados que forma la nueva malla, 
 !mayor, igual o menor a la anterior.
 rewind(2)
 do u=1,uu
 read(2,26)no(1),no(4),no(2),no(5),no(3),no(6)
   do uuu=1,6
   eval(no(uuu))=0.0_8
   enddo									    
 enddo

!Se aplica la condición seco-mojado (condición completa):
!--------------------------------------------------------
!Siempre necesaria cuando se hagan al menos 2 iteraciones consecutivas de la ecuación de aguas someras (siempre en el modelo superficial). 
if (modelo.eq.'superficial') then  
 !Tercer proceso. Se seleccionan otros elementos que también pertenecerían a la malla (hay calados positivos en los tres nodos esquina).
 rewind(4) 
 read(4,23)a
 do u=1,j
 read(4,26)no(1),no(4),no(2),no(5),no(3),no(6)
   if ((eval(no(1)).eq.0.0_8).and.(eval(no(2)).eq.0.0_8).and.(eval(no(3)).eq.0.0_8)) then	
	do ui=1,3
    b=ui+1-sb(ui)
    c=ui+2-sb(ui+1)
    e=ui+4-sb(ui)
    f=ui+5-sb(ui+1) 
    if (eval(no(ui+3)).ne.0.0_8) then
	uu=uu+1
	write(2,27)uu,no(1),no(4),no(2),no(5),no(3),no(6)
	!Se considera el elemento (sólo necesario dar valores a 'eval' en nodos cuadráticos).
	eval(no(ui+3))=0.0_8
	eval(no(e))=0.0_8
	eval(no(f))=0.0_8	
	 !Se darán dos casos bajo esta condición:
	 !if ((si(2*i+no(ui)).ne.0).and.(si(2*i+no(b)).ne.0).and.(si(2*i+no(c)).ne.0)) then 	 
	  !Selección de elementos encerrados.
	  !Al añadir elementos pegados al contornos móvil en zonas donde éste tiene gran curvatura pueden quedar atrapados 
	  !elementos de este tipo entre el nuevo contorno móvil. Se localizará un contorno (con el algoritmo donde se da vuc=0) que los tiene en cuenta. 	 
	  !Si no se seleccionasen se trabajaría con una malla de menos elementos donde algunos nodos del contorno no tendrán condición.
	  !Proceso necesario aunque supondrá un error si hay alturas muy distintas en los nodos de cada elemento a añadir.
	 !continue
	 !else
	  !Selección de otros elementos particulares.
	  !Pueden aparecer varios elementos con un nodo esquina común fuera del contorno móvil original de modo que sólo se añade alguno de ellos. 
	  !Proceso no necesario, pues siempre se dará condición en el contorno de la malla seleccionada.	   
	  !Se añaden elementos con poco sentido físico cuando la discretización no es buena (poco refinamiento). Siempre supondrá 
	  !un error ya que la altura será muy distinta en los nodos de cada elemento a añadir. Una posible solución en este caso es subir 'hmin'.
	 !continue 
	 !endif
	exit
	endif
	enddo
   endif									    
 enddo

 !Se sobreescribe 'vuc' para dar condición de contorno de no deslizamiento (velocidades nulas).  
 rewind(4)
 read(4,23)a
 do u=1,j
 read(4,26)no(1),no(4),no(2),no(5),no(3),no(6)				 
   !Aquí se localizan los nuevos elementos interfaz. Considerando que los calados pueden haber sido modificados, éstos serán  
   !Los elementos donde hay dos nodos de calado positivo y uno de calado negativo (uno o más elementos considerados en la iteración anterior ahora 
   !se eliminan). Único tipo de elementos si es la primera vez que se selecciona un dominio.
   !Los elementos donde haya dos nodos de calado positivos y uno de calado nulo (elemento ya considerado en la malla anterior o elemento añadido).
   !El objetivo es localizar los nodos del contorno móvil nuevo. Por ello la localización de elementos se hace mejor a través de 'eval'. 
   !Con ello se asegura que se están tratando nodos del contorno móvil nuevo. 
   do ui=1,3
   b=ui+1-sb(ui)
   c=ui+2-sb(ui+1)
   if ((eval(no(ui)).ne.sqrt(3.0_8)).and.(eval(no(b)).ne.sqrt(3.0_8)).and.(eval(no(c)).eq.sqrt(3.0_8))) then	 						 
	vuc(no(ui))=0.0_8   		  	 
    vuc(no(b))=0.0_8
	vuc(no(ui+3))=0.0_8
	vuc(i+no(ui))=0.0_8   		  	 
    vuc(i+no(b))=0.0_8
	vuc(i+no(ui+3))=0.0_8   	 	   
   endif
   enddo 
 enddo 

 !Cuarto proceso. Arreglo de condiciones de contorno en elementos casi despegados (pegados por un nodo) ya considerados. 
 !Es posible añadir un elemento (o más) tal que el nuevo contorno móvil haga una curva cerrada con él (o ellos). 
 !Esto puede ocurrir para un (dos o más) elemento si hay al menos tres (cuatro o más) elementos fuera del dominio original que compartan un nodo del 
 !contorno móvil y no tienen ningún otro nodo en este contorno. Si se tienen tres elementos así y se selecciona el elemento del medio (con los pasos
 !anteriores) aparecerá un problema. Con otra disposición geométrica el tercer proceso añadirá elementos a los lados de este elemento.  
 !El problema es que el elemento tendrá vuc=0 en todos sus nodos. Así, tanto el gradiente de la velocidad como la velocidad serán nulos y en la ecuación 
 !de continuidad no se aportan coeficientes en las filas y columnas correspondientes a los nuevos dos nodos (al menos sin estabilización).   
 do u=1,i
 pa(u)=0
 pe(u)=0
 enddo
 rewind(4)
 read(4,23)a
 !Selección de los nodos que causan el problema (para uno o más elementos de este tipo).
 do u=1,j
 read(4,26)no(1),no(4),no(2),no(5),no(3),no(6)				 
   if ((vuc(no(1)).eq.0.0_8).and.(vuc(no(2)).eq.0.0_8).and.(vuc(no(3)).eq.0.0_8)) then
   if ((vuc(no(4)).eq.0.0_8).and.(vuc(no(5)).eq.0.0_8).and.(vuc(no(6)).eq.0.0_8)) then
   !Si todos los nodos del elemento tienen condición de no delizamiento se sobrescribe la variable 'pa' en los nodos esquina.  
   !De momento la variable 'vuc' sólo es nula en los nodos del contorno móvil	 						 
	do ui=1,3
    pa(no(ui))=1
    enddo
	!Se evita que este elemento sea considerado como parte de la malla durante este proceso. Así, se trata una malla modificada.  	 	   
    cycle
   endif
   endif
   if ((eval(no(1)).eq.0.0_8).and.(eval(no(2)).eq.0.0_8).and.(eval(no(3)).eq.0.0_8)) then
   !Si el elemento está en la malla modificada se sobreescribe la variable 'pe' en los nodos esquina.
   do ui=1,3
   pe(no(ui))=1
   enddo
   endif 
 enddo
 do u=1,i
 !El problema se resuelve dando condiciones de altura en los dos nuevos nodos citados. Así no existirán esas ecuaciones elementales de continuidad.
 if ((pa(u).eq.1).and.(pe(u).eq.0)) then
 !Aquí se seleccionan los nuevos nodos. Para entrar aquí cada nodo debe pertenecer a un elemento con condición de contorno de no deslizamiento 
 !en todos los nodos y debe pertenece a un elemento que no está en la malla modificada.
 !Un valor nulo de condición de altura en vez de el vb(2*i+u), supuesto como valor de la iteración anterior, puede dar problemas de deselección 
 !y selección de estos elementos. 
 vuc(2*i+u)=vb(2*i+u) 		  	 
 endif
 enddo
 
 !En el nuevo contorno móvil se dan valores nulos de velocidad. Así, se calcularán bien las integrales de contorno.
 do u=1,2*i
  if (vuc(u).ne.sqrt(2.0_8)) then		   
  vt(u)=vuc(u)
  endif
 enddo

 !Fuera de la nueva malla considerada (en todo el dominio que no es de aguas someras) se dan condiciones de contorno nulas 
 !de velocidad y altura (u, v, Ht). Así, se reducirá el sistema eliminando filas y columnas que ya serán nulas (no se generarán coeficientes para
 !elementos que no son de esta malla). Finalmente no se calculará la solución en esos nodos.   
 do u=1,i
  if (eval(u).eq.sqrt(3.0_8)) then
  vuc(u)=0.0_8
  vuc(i+u)=0.0_8
  vuc(2*i+u)=0.0_8
  endif
 enddo

 !Introducción de las CC de flujo superficial que existen en el contorno no móvil (leídas desde fichero).
 do u=1,i+2
  read(4,'(A)')a
 enddo
 eaf=0
 do while (eaf.eq.0)
  read (4,46,iostat=eaf) u,fe,fi
  io=io+1  
  if (eaf.eq.0) then
   if (io.le.i) then
    if ((vuc(u).eq.sqrt(2.0_8)).and.(fe.ne.'o')) then
    vuc(u)=fi	   
    endif
   elseif ((io.gt.i).and.(io.le.2*i)) then
    if ((vuc(i+u).eq.sqrt(2.0_8)).and.(fe.ne.'o')) then
    vuc(i+u)=fi
    endif															   
   else
   !Al no darse condiciones de contorno de calado (vuc(2*i+u).eq.sqrt(2.0_8)) es equivalente a (eval(u).ne.sqrt(3.0)).
    if ((vuc(2*i+u).eq.sqrt(2.0_8)).and.(fe.ne.'o')) then	
    vuc(2*i+u)=fi+z(u)
    endif
   endif
  endif 
 enddo

 !Escritura de la nueva malla:  
 rewind(2)                             
 close(4)           
 open(1,file = 'C:\malla.txt',status='replace')
 write(1,'(A)')'nudos'
 write(1,10)i
 write(1,'(A)')'elementos'
 write(1,10)uu
 write(1,'(A)')'elemento p1  p2  p3  p4  p5  p6 (antihorarios)'
 do u=1,uu
  read(2,'(A)')a
  write(1,'(A)')a
 enddo
 close(2)
 close(1)
 j=uu

!Se aplica la condición similar a la condición seco mojado:
!----------------------------------------------------------
!Se aplica en otro caso (modelo conjunto). Sólo considera el primer y el segundo procesos.
!Con la subrutina mallaaguassomeras se llega por otro proceso a la selección que produce el tercer proceso de esta subrutina.
!No es necesario arreglar 'vuc' ya que se van a calcular velocidades subterráneas. Éstas se introducen en la subrutina mallaaguassomeras.
!Al utilizar la subrutina mallaaguassomeras habrá velocidades subterráneas en todos los nodos por lo que no será necesario el cuarto proceso
!de esta subrutina.
!Tampoco será necesario arreglar las condiciones de contorno o escribir la malla ya que se hará en la subrutina mallaaguassomeras.
else
 close(2)                             
 close(4)
 !La subrutina mallasubterranea se basará para localizar el contorno móvil (el mismo que localizaría la condición seco-mojado) en los 
 !valores nulos o negativos de calado. Podría ocurrir que quedasen nodos con calado positivo fuera del nuevo contorno móvil. 
 !Por ejemplo si en el primer proceso de esta subrutina no se considera un elemento que tiene dos nodos con calado positivo y uno de ellos en el 
 !nuevo contorno móvil. Por tanto:
 do u=1,i
  if ((eval(u).eq.sqrt(3.0_8)).and.(vt(2*i+u).gt.0.0_8)) then
  vt(2*i+u)=0.0_8
  endif
 enddo
 !Obviamente, esto no es necesario si se calcula a continuación otra iteración de la ecuación de aguas someras (con la condición seco-mojado).
endif
     
deallocate(si,pa,pe)
end

!---------------------------------------------------------------------------------------------------------------------------------------------
!Subrutina MALLASUBTERRANEA
!Esta subrutina puede modificar la malla del dominio subterráneo (sobreescribiendo el fichero mallasub.txt) haciendo una selección del dominio 
!subterráneo en base al dominio superficial. 
!Además da CC de altura, esto es, de nivel freático en el contorno móvil.
!Se aplica antes de resolver la ec subterránea y es equivalente a la subrutina mallaaguassomeras.
!Sólo se utiliza si se resuelve el modelo conjunto.
!Se utiliza una vez por cada grupo de iteraciones hasta convergencia del modelo subterráneo (en cada iteración conjunta) si se aplica el 
!modelo cunjunto. 
!---------------------------------------------------------------------------------------------------------------------------------------------
subroutine mallasubterranea (i,j,z,zp,zzp,vt,eval,qx,qy,vuc)
use interaccion
integer*4 i,j,u,ui,uuu,uu,io		   
real*8 vt(3*i),z(i),zp(i),zzp(i),eval(i),vuc(i),qx(i),qy(i),fi
character fe*1											 

10  format(I5)
21  format(3/,I5)								
23  format(4/,A80)									
26  format(6X,6(X,I5))	
27  format(I5,X,6(X,I5))	
46  format(X,I5,X,A1,X,F11.7)						
 
!Inizialización previa de variables:
!-----------------------------------
io=0
uu=0						  
do u=1,i					 
qx(u)=sqrt(3.0_8)					 
qy(u)=sqrt(3.0_8)
eval(u)=sqrt(3.0_8)
vuc(u)=sqrt(2.0_8)
zzp(u)=zp(u)
enddo

!Selección de dominio considerando los valores de calado:
!--------------------------------------------------------
 open(unit=7,file='C:\mallasubinicial.txt',status='old')
 open(2,status= 'scratch') 
 read(7,21)j
 read(7,'(A)')a
 do u=1,j
 read(7,26) no(1),no(4),no(2),no(5),no(3),no(6)
   !Selección inicial. Aquí sólo se consideran los elementos con todos sus nodos de calado nulo.  
   if ((vt(2*i+no(1)).le.0.0).and.(vt(2*i+no(2)).le.0.0).and.(vt(2*i+no(3)).le.0.0)) then
   uu=uu+1
   write(2,27)uu,no(1),no(4),no(2),no(5),no(3),no(6)
   do uuu=1,6
   eval(no(uuu))=0.0_8
   enddo
   endif
 enddo
 rewind(7)
 read(7,23)a
 do u=1,j
 read(7,26)no(1),no(4),no(2),no(5),no(3),no(6)
   do ui=1,3
   b=ui+1-sb(ui)
   c=ui+2-sb(ui+1)      
   !Se toman los elementos que rodean a la selección inicial (en ellos no se ha calculado el flujo superficial). También se toman elementos pegados 
   !a nodos aislados o a líneas que no forman elementos (regueros) donde no se calculó solución superficial (calados nulos). Es probable 
   !encontrar regueros pegados a la selección inicial. Se trata de los elementos interfaz (definidos en la subrutina nuevamalla).
   if ((vt(2*i+no(ui)).gt.0.0_8).and.(vt(2*i+no(b)).le.0.0_8).and.(vt(2*i+no(c)).le.0.0_8)) then  
	!Elementos con un nodo apoyado en el contorno móvil - tienen nodos esquina con: un calado positivo y dos nulos o negativos.
	uu=uu+1
    write(2,27)uu,no(1),no(4),no(2),no(5),no(3),no(6)
	do uuu=1,6
    eval(no(uuu))=0.0_8
    enddo
   elseif ((vt(2*i+no(ui)).gt.0.0_8).and.(vt(2*i+no(b)).gt.0.0_8).and.(vt(2*i+no(c)).le.0.0_8)) then
   !Elementos con dos nodos apoyados en el contorno móvil - tienen nodos esquina con: dos calado positivos y uno nulo o negativo.
   !No hay que considerar elementos con dos calados nulos y uno positivo para dar CC pues aún en caso de que haya dos o más elementos juntos de 
   !este tipo compartirán este nodo de calado positivo, y este nodo de la malla formará finalmente parte de un elemento de los analizados aquí. 
	uu=uu+1
    write(2,27)uu,no(1),no(4),no(2),no(5),no(3),no(6)
	do uuu=1,6
    eval(no(uuu))=0.0_8
    enddo
	!Se llega a una interfaz sumando elementos a la selección inicial. Será la interfaz que se calcule en la subrutina mallaaguassomeras.
	!Por tanto la interfaz será la misma para los dominios superficial y subterráneo. Por tanto se pasan valores que ya tienen esos nodos (dados
	!en la subrutina nuevamalla). Serán CC de altura.
	vuc(no(ui))=vt(2*i+no(ui))+z(no(ui))
	vuc(no(b))=vt(2*i+no(b))+z(no(b))   		  	 
	!Se modifica la cota del sustrato impermeable para garantizar que el flujo que entra o sale por el contorno móvil para el modelo subterráneo
	!sea el que sale o entra por el contorno móvil para el modelo superficial. 
	zzp(no(ui))=z(no(ui))
	zzp(no(b))=z(no(b))	      
   endif
   enddo   
 enddo

!Arreglo de las condiciones de contorno:
!---------------------------------------
!Fuera de la nueva malla considerada (en todo el dominio que no es subterráneo) se dan condiciones de contorno nulas 
!de nivel freático (hd). Así no se calculará la solución en estos nodos.	
do u=1,i
 if (eval(u).eq.sqrt(3.0_8)) then	   
 vuc(u)=0.0_8
 endif
enddo
	
!Introducción de las CC de flujo subterráneo que existen en el contorno no móvil (leídas desde fichero).    
do u=1,i+2
 read(7,'(A)')a
enddo
eaf=0
do while (eaf.eq.0)
 read (7,46,iostat=eaf) u,fe,fi
 io=io+1
 
!En la intersección del contorno móvil con el contorno del dominio conjunto también se tienen las CC de flujo subterráneo para el contorno del 
!dominio conjunto. Los nodos de dicha intersección pertenecerán a la zona de CC de flujo superficial del dominio conjunto (el contorno móvil está 
!en la orilla del dominio superficial). Así, en el fichero para flujo subterráneas se tiene CC de velocidad nula en ellos. 
!Por este motivo sólo se utiliza el valor dado de altura de la lámina en la intersección.
if (eaf.eq.0) then  
  if (io.le.i) then
	if ((vuc(u).eq.sqrt(2.0_8)).and.(fe.ne.'o')) then    
	qx(u)=fi
	endif
  elseif ((io.gt.i).and.(io.le.2*i)) then	 
	if ((vuc(u).eq.sqrt(2.0_8)).and.(fe.ne.'o')) then    
	qy(u)=fi
	endif															   
  else
    !En el contorno móvil, vuc será distinto de sqrt(2.0) ya que se impuso CC. Así, no se seleccionan los nodos de la intersección.
	if ((vuc(u).eq.sqrt(2.0_8)).and.(fe.ne.'o')) then 	   
	vuc(u)=fi
	endif
  endif
 endif 
enddo
		
!Escritura de la malla:
!----------------------
!Se escribe una nueva malla para subterráneo que no se va a modificar durante la resolución de la ecuación de agua subterránea.
!La malla tendrá la i original porque los números de nodo irán hasta i, y la nueva j para construir bien las cajas existentes.
rewind(2)                             
close(7)           
open(3,file = 'C:\mallasub.txt',status='replace')
write(3,'(A)')'nudos'
write(3,10)i
write(3,'(A)')'elementos'
write(3,10)uu
write(3,'(A)')'elemento p1  p2  p3  p4  p5  p6 (antihorarios)'
do u=1,uu
 read(2,'(A)')a
 write(3,'(A)')a
enddo
close(2)
close(3)
j=uu
end

!----------------------------------------------------------------------------------------------------------------------------------
!Subrutinas para el cálculo de las dimensiones de los vectores que llevan la matriz del sistema (predimensionamiento).
!----------------------------------------------------------------------------------------------------------------------------------
!Subrutina DIMVECTAS
!Cálculo de la dimensión para el sistema de ecuaciones para flujo superficial.
!Se consideran 9 posiciones en la matriz del sistema. La primera posición iría desde la fila 1 a la 'i' y de la columna 1 a la 'i'. 
!----------------------------------------------------------------------------------------------------------------------------------		  
subroutine dimvectas (i,j,newton,sino,vn,est)
use elemental
use allocatacion
integer*4, dimension(:),allocatable::nt,Eas,Easdosi,Eastresi,Eeqe,Eqqe,Eqe,Eee
integer*4 i,j,vn(i),Et
character newton*2,sino*2,est*2 

allocate(nt(i),Eas(3*i),Easdosi(3*i),Eastresi(3*i),Eeqe(i),Eqqe(i),Eqe(i),Eee(i))

39     format(4/,A80)
40     format(6X,6(X,I5))

!Inizialización previa de variables:
!-----------------------------------
Et=0
do u=1,i
nt(u)=0
Eqqe(u)=0
Eqe(u)=0
Eeqe(u)=0
Eee(u)=0
enddo

!Cálculo de la dimensión necesaria sin ensamblar para una caja determinada
!-------------------------------------------------------------------------
open(unit=1,file='C:\malla.txt',status='old')
read(1,39)ac
do u=1,j
read(1,40) n(1),n(4),n(2),n(5),n(3),n(6) 
 !Cálculo en 'nt(i)' del número de elementos que tiene el nodo 'i' en común.
 do ui=1,6
 nt(n(ui))=nt(n(ui))+1
 enddo
enddo
close(1)

!Se calcula el número de coeficientes por fila que se podrían generar. 
!'Eeqe' será el nº de coeficientes en los nodos esquina si se evalúa una función en nodos esquina y cuadráticos, 'Eqqe' el nº en los 
!nodos cuadráticos si se evalúa en nodos esquina y cuadráticos, 'Eee' el nº en los nodos esquina si se evalúa en los nodos esquina.
!'Eqe' el nº en los nodos cuadráticos si se evalúa en los nodos esquina.   
do u=1,i
 if (vn(u).eq.0) then	
 !Si 'u' es nodo esquina. 'Eqe(u)' y 'Eqqe(u)' serán nulos. 
  Eeqe(u)=6*nt(u)
  Eee(u)=3*nt(u) 
 else
 !Si 'u' no es nodo esquina. 'Eee(u)' y 'Eeqe(u)' serán nulos.					 
  Eqqe(u)=6*nt(u)
  Eqe(u)=3*nt(u)
 endif
enddo
!Así se pueden conocer el nº de coeficientes por fila para una caja determinada. Sólo se usaría Eeqe si se tuviese 
!una caja con funciones peso discretizadas en elemento lineal y función a integrar en elemento cuadrático.

!Cálculo de la dimensión para la matriz del sistema suponiendo que sólo se añada una caja allí donde haya cajas no nulas 
!-----------------------------------------------------------------------------------------------------------------------
!Sólo sería necesario calcular la variable 'Eas'.
if ((newton.eq.'no').and.(est.eq.'no')) then
!En este caso en las posiciones de la matriz del sistema (1,2) y (2,1) no se calcula ninguna caja. 
 do u=1,i
 !Se tienen en cuenta los tipos de cajas que hay y se calcula el nº total de coeficientes por fila.
 !En 'Eas(u)' está el nº de coeficientes por fila 'u'. Desde la fila 1 a la 'i':
 Eas(u)=Eeqe(u)+Eqqe(u)+Eee(u)+Eqe(u)
 !Desde la fila 'i+1' a la '2*i':
 Eas(u+i)=Eeqe(u)+Eqqe(u)+Eee(u)+Eqe(u)
  !Desde la fila '2*i+1' a la '3*i': 
  if (sino.eq.'si')then
  Eas(u+2*i)=2*Eeqe(u)+Eee(u)
  else
  !No se calcula caja para la posición (3,3).
  Eas(u+2*i)=2*Eeqe(u)
  endif
 !Se tienen en cuenta los tipos de cajas que hay y se calcula el nº de coeficientes hasta la columna 'i' por fila.
 Easdosi(u)=Eeqe(u)+Eqqe(u)
 Easdosi(u+i)=0
 Easdosi(u+2*i)=Eeqe(u)
 !Se tienen en cuenta los tipos de cajas que hay y se calcula el nº de coeficientes hasta la columna '2*i' por fila.
 Eastresi(u)=Eeqe(u)+Eqqe(u)
 Eastresi(u+i)=Eeqe(u)+Eqqe(u)
 Eastresi(u+2*i)=2*Eeqe(u)
 enddo
else	
!En otro caso se considera que hay cajas en todas las posiciones.
 do u=1,i
 !Se calcula el nº total de coeficientes por fila.
 Eas(u)=2*Eeqe(u)+2*Eqqe(u)+Eee(u)+Eqe(u)
 Eas(u+i)=2*Eeqe(u)+2*Eqqe(u)+Eee(u)+Eqe(u)
 Eas(u+2*i)=2*Eeqe(u)+Eee(u)
 !Se calcula el nº de coeficientes hasta la columna 'i' por fila.
 Easdosi(u)=Eeqe(u)+Eqqe(u)
 Easdosi(u+i)=Eeqe(u)+Eqqe(u)
 Easdosi(u+2*i)=Eeqe(u)
 !Se calcula el nº de coeficientes hasta la columna '2*i' por fila.
 Eastresi(u)=2*Eeqe(u)+2*Eqqe(u)
 Eastresi(u+i)=2*Eeqe(u)+2*Eqqe(u)
 Eastresi(u+2*i)=2*Eeqe(u)
 enddo
endif
!Cálculo de la dimensión total que tendrán los vectores. Es necesario reservar un espacio de '3*i+1' para referenciar 
!en que posición empiezan los coeficientes de cada fila (se hace a continuación).
do u=1,i
Et=Et+Eas(u)+Eas(u+i)+Eas(u+2*i)
enddo
Et=Et+3*i+1

!Se referencia el vector puntero
!-------------------------------
allocate (cia(Et),ca(Et),ita(Et),sa(Et))
do u=1,Et
ita(u)=0
sa(u)=0.0_8
enddo
!Se referencia en cada componente 'u' del vector en que posición se empezarán a escribir los coeficientes de la fila 'u'.
!La suma del nº total de coeficientes por fila 'Eas' dará información de la posición del primer coeficiente de la siguiente fila.  
!Finalmente 'ita(3*i+1)'-1 es igual a 'Et' (posteriormente ndim).
ita(1)=3*i+2	   
do u=1,3*i
ita(u+1)=ita(u)+Eas(u)
enddo

!Cálculo de otras variables para no reservar más coeficientes para cajas que se sumen sobre otras.
!-------------------------------------------------------------------------------------------------
!Se calcula la variable 'posi(u)' con posición ocupada en 'ita' por el primer coeficiente de cada fila 'u'. También la variable 
!'posdosi(u)' con posición ocupada en 'ita' por el primer coeficiente situado en una columna posterior a 'i' de cada fila 'u'. 
!La suma del nº de coeficientes por fila 'Easdosi' dará esta información. También la variable 'postresi(u)' con posición ocupada 
!en 'ita' por el primer coeficiente situado en una columna posterior a '2*i' de cada fila 'u'.
do u=1,3*i
posi(u)=ita(u)
posdosi(u)=ita(u)+Easdosi(u)
postresi(u)=ita(u)+Eastresi(u)
enddo
!Al escribir cajas en mismos lugares basta con inicializar posiciones utilizando estas variables y sumar los nuevos coeficientes 
!a los existentes. Esto se puede hacer dado que se crean coeficientes cuya referencia sigue el mismo orden (los elementos se leen 
!por el mismo orden). 

deallocate(nt,Eas,Easdosi,Eastresi,Eeqe,Eqqe,Eqe,Eee)
end

!---------------------------------------------------------------------------------------------------------------------------------
!Subrutina DIMVECTSB
!Cálculo de la dimensión para el sistema de ecuaciones para flujo subterráneo.
!Se considera que solo hay 1 posición en la matriz del sistema. 
!---------------------------------------------------------------------------------------------------------------------------------
subroutine dimvectsb (i,j,vn)
use elemental
use allocatacion
integer*4, dimension(:),allocatable::nt,Esb,Eee
integer*4 i,j,vn(i),Et	

allocate(nt(i),Esb(i),Eee(i))

39     format(4/,A80)
40     format(6X,6(X,I5))

!Inizialización previa de variables:
!-----------------------------------
Et=0
do u=1,i
nt(u)=0
Eee(u)=0
enddo

!Cálculo de la dimensión necesaria sin ensamblar para una caja determinada
!-------------------------------------------------------------------------
open(unit=1,file='C:\mallasub.txt',status='old')
read(1,39)ac
do u=1,j
read(1,40) n(1),n(4),n(2),n(5),n(3),n(6) 
 do ui=1,6
 nt(n(ui))=nt(n(ui))+1
 enddo
enddo
close(1)
    
do u=1,i
 if (vn(u).eq.0) then
 !Si 'u' es nodo esquina. Eee(u) será nulo si 'u' no es esquina.	 
  Eee(u)=3*nt(u) 
 endif
enddo
!Sólo habrá caja con funciones peso discretizadas en elemento lineal y función a integrar en elemento lineal.

!Cálculo de la dimensión para la matriz del sistema suponiendo que sólo se añada una caja 
!----------------------------------------------------------------------------------------
do u=1,i
Esb(u)=Eee(u)
enddo
do u=1,i
Et=Et+Esb(u)
enddo
Et=Et+i+1

!Se referencia el vector puntero
!-------------------------------
allocate (cia(Et),ca(Et),ita(Et),sa(Et))
do u=1,Et
ita(u)=0
sa(u)=0.0_8
enddo
ita(1)=i+2	   
do u=1,i
ita(u+1)=ita(u)+Esb(u)
enddo

!Cálculo de otras variables para no reservar más coeficientes para cajas que se sumen sobre otras.
!-------------------------------------------------------------------------------------------------
!Se calcula una variable 'pos(u)' con posición ocupada en 'ita' por el primer coeficiente de cada fila 'u'.
do u=1,i
pos(u)=ita(u)
enddo

deallocate(nt,Esb,Eee)
end

!----------------------------------------------------------------------------------------------------------------------------------
!Subrutinas para el cálculo de coeficientes de las integrales de contorno (vectores elementales).
!----------------------------------------------------------------------------------------------------------------------------------
!Subrutina VECTORCONTORNOPRESIONES
!En esta subrutina se calculan las integrales de contorno para la ecuación de Navier-Stokes 2D o para la ecuación de aguas someras.
!La integración se ha hecho analíticamente (a mano) obteniendo las expresiones necesarias. 
!----------------------------------------------------------------------------------------------------------------------------------
subroutine vectorcontornopresiones(i,j,x,y,vv,vvv,vectorx,vectory,vector,nu,ten)
integer*4 i,j,u,ui,sb(4),n(6),a,b,c,d,e,f
real*8 vectorx(i),vectory(i),vector(i),vv(3*i),vvv(3*i),x(i),y(i),xa,xb,ya,yb,pu,pv,qu,qv,nu,ax,ay,ma,mb,mc,md,jac            
character ac*80,ten*2
 
39     format(4/,A80)
40     format(6X,6(X,I5))

!Inizialización previa de variables:
!-----------------------------------
ma=3.0_8/20.0_8
mb=1.0_8/60.0_8
mc=1.0_8/5.0_8
md=2.0_8/15.0_8
sb=(/0,0,3,3/)

!Cálculo de los valores para cada lado del elemento
!-------------------------------------------------- 
!Ensamblaje directo al sobreescribir los vectores con la suma de su valor y el de los coeficientes generados.
open(unit=1,file='C:\malla.txt',status='old')
read(1,39)ac
do u=1,j
!Selección de cada elemento
read(1,40) n(1),n(4),n(2),n(5),n(3),n(6)
jac=abs((x(n(1))-x(n(3)))*(y(n(2))-y(n(3)))-(x(n(2))-x(n(3)))*(y(n(1))-y(n(3))))
 do ui=1,3
  !Selección de cada lado del elemento.
  a=ui
  b=ui+1-sb(ui)
  c=ui+2-sb(ui+1)
  d=ui+3
  e=ui+4-sb(ui)
  f=ui+5-sb(ui+1)
 
 !Cálculo previo de coeficientes necesarios
 !-----------------------------------------
 !Coeficientes calculados por separado por la longitud de la expresión.     	
 pu=(vv(n(c))*ma-vv(n(b))*mb+vv(n(e))*mc)*vv(2*i+n(c))+(vv(n(c))*mb+vv(n(b))*mb+vv(n(e))*md)*vv(2*i+n(b))
 pv=(vv(i+n(c))*ma-vv(i+n(b))*mb+vv(i+n(e))*mc)*vv(2*i+n(c))+(vv(i+n(c))*mb+vv(i+n(b))*mb+vv(i+n(e))*md)*vv(2*i+n(b)) 
 qu=(vv(n(c))*mb+vv(n(b))*mb+vv(n(e))*md)*vv(2*i+n(c))+(-vv(n(c))*mb+vv(n(b))*ma+vv(n(e))*mc)*vv(2*i+n(b))
 qv=(vv(i+n(c))*mb+vv(i+n(b))*mb+vv(i+n(e))*md)*vv(2*i+n(c))+(-vv(i+n(c))*mb+vv(i+n(b))*ma+vv(i+n(e))*mc)*vv(2*i+n(b))
 !Cálculo de los catetos correspondientes al producto de la longitud del lado por coseno o seno del ángulo que forman la dirección normal 
 !y las direcciones de los ejes. 
 xa=x(n(c))-x(n(b)) 
 xb=x(n(a))-x(n(b))		   
 ya=y(n(c))-y(n(b))		
 yb=y(n(a))-y(n(b))
 !Cálculo del signo de las componentes del vector normal al contorno que tiene orientación saliente.
 call direccion(x,y,i,n(b),n(c),n(a),ax,ay)
 
 !Para integrales de contorno de las ecuaciones dinámicas - término de presión
 !----------------------------------------------------------------------------
 vectorx(n(c))=vectorx(n(c))-(vvv(2*i+n(c))/6.0_8)*ax*abs(ya)*9.81_8 
 vectorx(n(b))=vectorx(n(b))-(vvv(2*i+n(b))/6.0_8)*ax*abs(ya)*9.81_8
 vectorx(n(e))=vectorx(n(e))-((vvv(2*i+n(b))+vvv(2*i+n(c)))/3.0_8)*ax*abs(ya)*9.81_8 
 vectory(n(c))=vectory(n(c))-(vvv(2*i+n(c))/6.0_8)*ay*abs(xa)*9.81_8 
 vectory(n(b))=vectory(n(b))-(vvv(2*i+n(b))/6.0_8)*ay*abs(xa)*9.81_8 
 vectory(n(e))=vectory(n(e))-((vvv(2*i+n(b))+vvv(2*i+n(c)))/3.0_8)*ay*abs(xa)*9.81_8 
 
 !Para integrales de contorno de las ecuaciones dinámicas - término viscoso	(tensiones)
 !-------------------------------------------------------------------------------------
 !El usuario decide si considerar estos integrales o no. Tras ensamblar generarán valores no nulos en los nodos internos de la malla
 !dado que no existe continuidad de la derivada de las funciones de interpolación en la dirección normal (sólo la habría en la dirección del lado).
 !Así, estos términos pueden dar lugar a una mayor inestabilidad.
 if (ten.eq.'si') then
 vectorx(n(c))=vectorx(n(c))+(3.0_8*vv(n(c))+vv(n(b))-4.0_8*vv(n(e)))*(yb*ax*abs(ya)-xb*ay*abs(xa))*nu/(6.0_8*jac)
 vectorx(n(b))=vectorx(n(b))+(-vv(n(c))-3.0_8*vv(n(b))+4.0_8*vv(n(e)))*(yb*ax*abs(ya)-xb*ay*abs(xa))*nu/(6.0_8*jac)
 vectorx(n(e))=vectorx(n(e))+4.0_8*(vv(n(c))-vv(n(b)))*(yb*ax*abs(ya)-xb*ay*abs(xa))*nu/(6.0_8*jac)
 vectory(n(c))=vectory(n(c))+(3.0_8*vv(i+n(c))+vv(i+n(b))-4.0_8*vv(i+n(e)))*(yb*ax*abs(ya)-xb*ay*abs(xa))*nu/(6.0_8*jac)
 vectory(n(b))=vectory(n(b))+(-vv(i+n(c))-3.0_8*vv(i+n(b))+4.0_8*vv(i+n(e)))*(yb*ax*abs(ya)-xb*ay*abs(xa))*nu/(6.0_8*jac)
 vectory(n(e))=vectory(n(e))+4.0_8*(vv(i+n(c))-vv(i+n(b)))*(yb*ax*abs(ya)-xb*ay*abs(xa))*nu/(6.0_8*jac)
 endif

 !Para integrales de contorno de la ecuación de continuidad
 !---------------------------------------------------------
 !pu=q referente al nodo n(c) de integrar el calado por u. pv=q referente al nodo n(c) de integrar el calado por v
 !qu=q referente al nodo n(b) de integrar el calado por u. qv=q referente al nodo n(b) de integrar el calado por v
 vector(n(c))=vector(n(c))+(pu*ax*abs(ya)+pv*ay*abs(xa))
 vector(n(b))=vector(n(b))+(qu*ax*abs(ya)+qv*ay*abs(xa))  
	
 enddo
enddo
close(1)
end

!--------------------------------------------------------------------------------------------------
!Subrutina VECTORCONTORNOCAUDALES
!En esta subrutina se calculan las integrales de contorno para la ecuación para flujo subterraneo.
!La integración se ha hecho analíticamente (a mano) obteniendo las expresiones necesarias.
!--------------------------------------------------------------------------------------------------
subroutine vectorcontornocaudales(i,j,x,y,vector,qxx,qyy)
integer*4 i,j,sb(4),n(6),u,ui,a,b,c
real*8 vector(i),x(i),y(i),xa,ya,ax,ay,qxx(i),qyy(i)            
character ac*80

39     format(4/,A80)
40     format(6X,6(X,I5))
		  
!Inizialización previa de variables:
!-----------------------------------
sb=(/0,0,3,3/)

!Cálculo de los valores para cada lado del elemento
!--------------------------------------------------
open(unit=3,file='C:\mallasub.txt',status='old')
read(3,39)ac
do u=1,j
read(3,40)n(1),n(4),n(2),n(5),n(3),n(6)
 do ui=1,3
   a=ui
   b=ui+1-sb(ui)
   c=ui+2-sb(ui+1)
   xa=x(n(c))-x(n(b)) 
   ya=y(n(c))-y(n(b)) 
   call direccion(x,y,i,n(b),n(c),n(a),ax,ay)
   
   !Cálculo de las Integrales de contorno de caudal	(valor interpolado linealmente)	que aparecen en la ecuación. 
   vector(n(c))=vector(n(c))+((2.0_8*qxx(n(c))+qxx(n(b)))*ax*abs(ya)+(2.0_8*qyy(n(c))+qyy(n(b)))*ay*abs(xa))/6.0_8
   vector(n(b))=vector(n(b))+((2.0_8*qxx(n(b))+qxx(n(c)))*ax*abs(ya)+(2.0_8*qyy(n(b))+qyy(n(c)))*ay*abs(xa))/6.0_8
   !En el contorno móvil se podrían calcular unas integrales de contorno más complicadas en función de velx y vely pero da igual ya 
   !que los valores generados desaparecen al imponer condiciones de nivel freático.
 enddo
enddo
close(3)
end

!------------------------------------------------------------------------------------------------------------------------------------
!Subrutina DIRECCION
!Esta subrutina obliga a que la orientación del vector normal en cada lado del contorno elemental sea saliente respecto al elemento.
!Cada vez que se utiliza se analiza un lado del elemento. Así, se analiza el lado formado por los nodos 'k' y 'm'.
!'a' es el signo de la componente horizontal del vector, 'b' es el signo de la componente vertical del vector.
!------------------------------------------------------------------------------------------------------------------------------------
subroutine direccion(x,y,i,k,m,p,a,b)
integer*4 i,k,m,p
real*8 x(i),y(i),a,b

!La dirección del lado es distinta a la horizontal o a la vertical.
if (((x(k)-x(m)).ne.0.0).and.((y(k)-y(m)).ne.0.0)) then
  if ((y(p)-y(k)-((y(k)-y(m))/(x(k)-x(m)))*(x(p)-x(k))).lt.0.0) then
  b=1.0_8
	if ((x(p)-x(k)-((x(k)-x(m))/(y(k)-y(m)))*(y(p)-y(k))).lt.0.0) then
	a=1.0_8
	elseif ((x(p)-x(k)-((x(k)-x(m))/(y(k)-y(m)))*(y(p)-y(k))).gt.0.0) then
	a=-1.0_8
	endif
  elseif ((y(p)-y(k)-((y(k)-y(m))/(x(k)-x(m)))*(x(p)-x(k))).gt.0.0) then
  b=-1.0_8
    if ((x(p)-x(k)-((x(k)-x(m))/(y(k)-y(m)))*(y(p)-y(k))).lt.0.0) then
    a=1.0_8
	elseif ((x(p)-x(k)-((x(k)-x(m))/(y(k)-y(m)))*(y(p)-y(k))).gt.0.0) then
	a=-1.0_8
	endif
  endif
!La dirección del lado es vertical.
elseif ((x(k)-x(m)).eq.0.0) then
b=0.0_8 
  if ((x(k)-x(p)).gt.0.0) then
  a=1.0_8
  elseif ((x(k)-x(p)).lt.0.0) then
  a=-1.0_8
  endif 
!La dirección del lado es horizontal.
elseif ((y(k)-y(m)).eq.0.0) then
a=0.0_8 
  if ((y(p)-y(k)).lt.0.0) then
  b=1.0_8
  elseif ((y(p)-y(k)).gt.0.0) then
  b=-1.0_8
  endif
endif 				 
end

!-----------------------------------------------------------------------------------------------------------------------------------------
!Subrutinas para el cálculo de coeficientes de las integrales en el dominio (formulación Bubnov-Galerkin).
!Las funciones de forma están escrita en coordenadas locales. 
!Para elementos lineales 'Mp' serán las funciones de forma, 'Mpx' las derivadas en x de las funciones de forma, 'Mpy' las derivadas en y 
!de las funciones de forma. Para elementos cuadráticos 'Mi' serán las funciones de forma, 'Mix' serán las derivadas en x de las funciones 
!de forma, 'Miy' las derivadas en y de las funciones de forma.
!-----------------------------------------------------------------------------------------------------------------------------------------
!Subrutina VECTORCONTORNOFUENTE
!En esta subrutina se imponen las condiciones interiores de lluvia para flujo superficial calculando vectores elementales 
!que irán al término independiente del sistema.	
!Se hace una integración en el dominio de la intensidad de lluvia (tiene unidades m3/m2*s), cuyos valores nodales están guardados en 
!la variable 'ql'. Las integrales son calculadas numéricamente mediante cuadraturas.
!La función intensidad se interpola sólo a través de los nodos esquina (interpolación lineal). Es irrelevante como se interpole ya que 
!los valores nodales serán idénticos e iguales a la intensidad que define el usuario para todo el dominio. 
!'n' es una variable global, tiene dimensión 6 y lleva los números de los nodos de cada elemento leídos desde fichero. 
!-----------------------------------------------------------------------------------------------------------------------------------------
subroutine vectorcontornofuente(i,j,x,y,vectorz,ql)
use elemental
integer*4 i,j
real*8 x(i),y(i),Ip,Mp(3),vectorz(i),ql(i)            

39     format(4/,A80)
40     format(6X,6(X,I5))

!Inizialización previa de variables:
!-----------------------------------
!En 'luno,ldos,ja', que son variables globales, van los puntos y pesos de integración para la cuadratura.
!Ciertos valores son nulos dado que su dimensión es 7 y sólo se usarán 3 puntos de integración. 
luno=(/0.5_8, 0.0_8, 0.5_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8/)		   
ldos=(/0.5_8, 0.5_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8/)
ja=(/1.0_8/6.0_8, 1.0_8/6.0_8, 1.0_8/6.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8/)

!Cálculo de los valores para cada elemento
!-----------------------------------------
!Se calculan valores sólo para elementos del dominio, que están guardados en malla.txt
open(unit=1,file='C:\malla.txt',status='old')
read(1,39)ac
 do u=1,j
 !Los nodos en el fichero están escritos teniendo en cuenta el sentido antihorario. Las funciones de forma para elementos lineales 
 !están definidas de forma que están referidas a los nodos esquina en sentido antihorario. Por ello se almacenan los nodos en este orden.
 read(1,40) n(1),n(4),n(2),n(5),n(3),n(6)

	 !Se calculan antes estos valores ya que son contantes para cada elemento. Siempre se calcularía lo mismo si se evalúan estas expresiones 
	 !para cada punto de integración. Daría igual meterlo en el siguiente bucle aunque conllevaría mayor coste computacional.
	 xa=x(n(1))-x(n(3)) 
	 xb=x(n(2))-x(n(3)) 
	 ya=y(n(1))-y(n(3))
	 yb=y(n(2))-y(n(3))  
	 jac=abs(xa*yb-xb*ya)	

	 !Integración con tres puntos. Polinomios de grado superior a grado 2 quedarían integrados con error.
	 do uu=1,3	      
	 !Funciones necesarias.
	 Mp(1)=luno(uu)
	 Mp(2)=ldos(uu)
	 Mp(3)=1.0_8-luno(uu)-ldos(uu)	 
	 
	  !Se conocen valores de intensidad en los nodos esquina.
	  Ip=Mp(1)*ql(n(1))+Mp(2)*ql(n(2))+Mp(3)*ql(n(3))
	  !Ensamblaje directo de los vectores elementales.
	  do uuu=1,3
	  vectorz(n(uuu))=vectorz(n(uuu))+(Mp(uuu)*Ip)*ja(uu)*jac 	  
	  enddo	
	 enddo
 enddo
close(1)
end

!---------------------------------------------------------------------------------------------------------------------------------------
!Subrutina VECTORCONTORNOFUENTESUB
!En esta subrutina se imponen las condiciones interiores de lluvia para flujo subterráneo calculando vectores elementales 
!que irán al término independiente del sistema.	
!Se hace una integración en el dominio de la intensidad de lluvia (tiene unidades m3/m2*s), cuyos valores nodales están guardados en 
!la variable 'ql'. 
!La función intensidad se interpola sólo a través de los nodos esquina (interpolación lineal). 
!---------------------------------------------------------------------------------------------------------------------------------------
subroutine vectorcontornofuentesub(i,j,x,y,vectorz,ql)
use elemental
integer*4 i,j
real*8 x(i),y(i),Ip,Mp(3),vectorz(i),ql(i)            

39     format(4/,A80)
40     format(6X,6(X,I5))

!Inizialización previa de variables:
!-----------------------------------
luno=(/0.5_8, 0.0_8, 0.5_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8/)		   
ldos=(/0.5_8, 0.5_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8/)
ja=(/1.0_8/6.0_8, 1.0_8/6.0_8, 1.0_8/6.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8/)

!Cálculo de los valores para cada elemento
!-----------------------------------------
!Se calculan valores sólo para elementos del dominio, que están guardados en mallasub.txt
open(unit=3,file='C:\mallasub.txt',status='old')
read(3,39)ac
 do u=1,j
 read(3,40) n(1),n(4),n(2),n(5),n(3),n(6)

	 xa=x(n(1))-x(n(3)) 
	 xb=x(n(2))-x(n(3)) 
	 ya=y(n(1))-y(n(3))
	 yb=y(n(2))-y(n(3))  
	 jac=abs(xa*yb-xb*ya)	

	 do uu=1,3		      
	 
	 Mp(1)=luno(uu)
	 Mp(2)=ldos(uu)
	 Mp(3)=1.0_8-luno(uu)-ldos(uu)
	 
	  !Se conocen valores de intensidad en los nodos esquina.
	  Ip=Mp(1)*ql(n(1))+Mp(2)*ql(n(2))+Mp(3)*ql(n(3))
	  do uuu=1,3
	  vectorz(n(uuu))=vectorz(n(uuu))+(Mp(uuu)*Ip)*ja(uu)*jac 	  
	  enddo	
	 enddo
 enddo
close(3)
end

!-----------------------------------------------------------------------------------------------------------------------------------------
!Subrutina CAJASAB
!A continuación se calculan coeficientes de las matrices elementales que forman las caja A, Bx y By.
!'Nn' son las matrices elementales de la caja A, 'Nm' son las matrices elementales de la caja Bx y 'Nk' son las matrices elementales 
!de la caja By. Aparecen en las ecuaciones dinámicas para flujo superficial.
!'sa' es el vector que lleva los coeficientes de la matriz del sistema no nulos de fuera de la diagonal y 'ita' es el vector puntero 
!que lleva su posición (aún sin formato MSR). Son variables globales que ya han sido definidas.
!'poi,podosi,potresi' son inicializadas con las variables globales 'posi,posdosi,postresi' para indicar en que posición del vector para 
!cada fila se empiezan a escribir los coeficientes. De esta forma, al escribir los coeficientes uno tras otro, se suman de forma 
!adecuada si se escriben dos cajas de matrices en la misma posición de la matriz del sistema.
!Por ejemplo 'poi' lo indica para las posiciones (1,1),(2,1) y (3,1); 'podosi' lo indica para las posiciones (1,2),(2,2) y (3,2);... 
!Los coeficientes de las matrices elementales se integran y se escriben sobre los vectores ita y sa a través de la subrutina suma.
!----------------------------------------------------------------------------------------------------------------------------------------- 
subroutine cajasab(i,j,x,y,nu,del)
use elemental
integer*4, dimension(:),allocatable::poi,podosi,potresi
integer*4 i,j
real*8 x(i),y(i),Nn(6,6),Nm(6,3),Nk(6,3),Mix(6),Miy(6),Mp(3),nu,del

allocate(poi(3*i),podosi(3*i),potresi(3*i))
   
39     format(4/,A80)
40     format(6X,6(X,I5))

!Inizialización previa de variables:
!-----------------------------------
luno=(/0.5_8, 0.0_8, 0.5_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8/)		   
ldos=(/0.5_8, 0.5_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8/)
ja=(/1.0_8/6.0_8, 1.0_8/6.0_8, 1.0_8/6.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8/)
do u=1,3*i
poi(u)=posi(u)
podosi(u)=posdosi(u)
potresi(u)=postresi(u)
enddo

!Cálculo de los valores para cada elemento
!-----------------------------------------
open(unit=1,file='C:\malla.txt',status='old')
read(1,39)ac
 do u=1,j
     !Los nodos en el fichero están escritos teniendo en cuenta el sentido antihorario.
	 !Las funciones de forma para elementos cuadráticos están definidas de forma que están referidas primero a los nodos esquina en sentido 
	 !antihorario y después a los tres nodos cuadráticos en sentido antihorario. Por ello se almacenan los nodos en este orden.	 
     read(1,40) n(1),n(4),n(2),n(5),n(3),n(6)
  
     xa=x(n(1))-x(n(3)) 
	 xb=x(n(2))-x(n(3)) 
	 ya=y(n(1))-y(n(3))
	 yb=y(n(2))-y(n(3))
	 jac=abs(xa*yb-xb*ya)

	 !Inicialización de las matrices elementales para el cálculo de nuevas matrices elementales.
	 do uu=1,6
	 do uuu=1,6
     Nn(uu,uuu)=0.0_8
	 enddo
	 do uuu=1,3
	 Nm(uu,uuu)=0.0_8
     Nk(uu,uuu)=0.0_8
     enddo
	 enddo

	 do uu=1,3
	 !Los valores generados en un punto de integración se suman a los generados previamente. Así serán calculadas las integrales 
	 !y los valores irán en las matrices elementales al salir de este bucle.
	 
	 Mix(1)=((4.0_8*luno(uu)-1.0_8)*yb)/jac
	 Mix(2)=((1.0_8-4.0_8*ldos(uu))*ya)/jac
	 Mix(3)=((4.0_8*luno(uu)+4.0_8*ldos(uu)-3.0_8)*(yb-ya))/jac
	 Mix(4)=(4.0_8*ldos(uu)*yb-4.0_8*luno(uu)*ya)/jac
	 Mix(5)=((8.0_8*ldos(uu)+4.0_8*luno(uu)-4.0_8)*ya-4.0_8*ldos(uu)*yb)/jac
	 Mix(6)=((4.0_8-4.0_8*ldos(uu)-8.0_8*luno(uu))*yb+4.0_8*luno(uu)*ya)/jac
	 
	 Miy(1)=((1.0_8-4.0_8*luno(uu))*xb)/jac
	 Miy(2)=((4.0_8*ldos(uu)-1.0_8)*xa)/jac
	 Miy(3)=((4.0_8*luno(uu)+4.0_8*ldos(uu)-3.0_8)*(xa-xb))/jac
	 Miy(4)=(4.0_8*luno(uu)*xa-4.0_8*ldos(uu)*xb)/jac
	 Miy(5)=(4.0_8*ldos(uu)*xb+(4.0_8-4.0_8*luno(uu)-8.0_8*ldos(uu))*xa)/jac
	 Miy(6)=((8.0_8*luno(uu)+4.0_8*ldos(uu)-4.0_8)*xb-4.0_8*luno(uu)*xa)/jac
	 
	 Mp(1)=luno(uu)
	 Mp(2)=ldos(uu)
	 Mp(3)=1.0_8-luno(uu)-ldos(uu)
	 		
	  do ui=1,6
  	  do uj=1,6
	  !Posicion 1,1 y 2,2. Caja A.
	  Nn(uj,ui)=Nn(uj,ui)+(Mix(uj)*Mix(ui)+Miy(uj)*Miy(ui))*ja(uu)*jac
	  enddo
	  enddo
	  do ui=1,3
  	  do uj=1,6   	 
	  !Posicion 1,3. Caja Bx. 
	  Nm(uj,ui)=Nm(uj,ui)+(Mix(uj)*Mp(ui))*ja(uu)*jac	 				   
	  !Posicion 2,3. Caja By.
	  Nk(uj,ui)=Nk(uj,ui)+(Miy(uj)*Mp(ui))*ja(uu)*jac
	  enddo
	  enddo  
    enddo
	!Escritura de los coeficientes de las matrices elementales en los vectores 'ita,sa' llamando a la subrutina suma.
	!La caja A se escribe en dos posiciones diferentes.
	call suma(nu*Nn/del,n,0,0,3*i,6,6,poi)
	call suma(nu*Nn/del,n,i,i,3*i,6,6,podosi)
	call suma(-Nm*9.81_8/del,n,0,2*i,3*i,6,3,potresi) 
	call suma(-Nk*9.81_8/del,n,i,2*i,3*i,6,3,potresi)
enddo
close(1)
deallocate(poi,podosi,potresi)
end

!----------------------------------------------------------------------------------------------------------------------------------
!Subrutina CAJASBT
!A continuación se calculan coeficientes de las matrices elementales que forman las cajas Bxt y Byt.
!'Nn' son las matrices elementales de la caja Bxt, 'Nm' son las matrices elementales de la caja Byt. Son las transpuestas de las
!cajas Bx y By definidas en la subrutina cajasab.
!Estas cajas aparecen al discretizar la ecuación de continuidad para las ecuaciones Navier-Stokes 2D. Estas ecuaciones pueden ser 
!utilizadas en la primera o primeras iteraciones de Picard para buscar una buena aproximación inicial para las ecuaciones de aguas 
!someras. 
!----------------------------------------------------------------------------------------------------------------------------------
subroutine cajasbt (i,j,x,y,del)
use elemental
integer*4, dimension(:),allocatable::poi,podosi
integer*4 i,j
real*8 x(i),y(i),Nn(3,6),Nm(3,6),Mix(6),Miy(6),Mp(3),del

allocate(poi(3*i),podosi(3*i))

 39     format(4/,A80)
 40     format(6X,6(X,I5))

!Inizialización previa de variables:
!-----------------------------------
luno=(/0.5_8, 0.0_8, 0.5_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8/)		   
ldos=(/0.5_8, 0.5_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8/)
ja=(/1.0_8/6.0_8, 1.0_8/6.0_8, 1.0_8/6.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8/)
do u=1,3*i
poi(u)=posi(u)
podosi(u)=posdosi(u)
enddo

!Cálculo de los valores para cada elemento
!-----------------------------------------
open(unit=1,file='C:\malla.txt',status='old')
read(1,39)ac
 do u=1,j 
     read(1,40) n(1),n(4),n(2),n(5),n(3),n(6)
 	 
	 xa=x(n(1))-x(n(3)) 
	 xb=x(n(2))-x(n(3)) 
	 ya=y(n(1))-y(n(3))
	 yb=y(n(2))-y(n(3))
	 jac=abs(xa*yb-xb*ya)
	   
	 do uu=1,3
     do uuu=1,6
     Nn(uu,uuu)=0.0_8
     Nm(uu,uuu)=0.0_8
     enddo
     enddo

	 do uu=1,3     
	 
	 Mix(1)=((4.0_8*luno(uu)-1.0_8)*yb)/jac
	 Mix(2)=((1.0_8-4.0_8*ldos(uu))*ya)/jac
	 Mix(3)=((4.0_8*luno(uu)+4.0_8*ldos(uu)-3.0_8)*(yb-ya))/jac
	 Mix(4)=(4.0_8*ldos(uu)*yb-4.0_8*luno(uu)*ya)/jac
	 Mix(5)=((8.0_8*ldos(uu)+4.0_8*luno(uu)-4.0_8)*ya-4.0_8*ldos(uu)*yb)/jac
	 Mix(6)=((4.0_8-4.0_8*ldos(uu)-8.0_8*luno(uu))*yb+4.0_8*luno(uu)*ya)/jac
	
	 Miy(1)=((1.0_8-4.0_8*luno(uu))*xb)/jac
	 Miy(2)=((4.0_8*ldos(uu)-1.0_8)*xa)/jac
	 Miy(3)=((4.0_8*luno(uu)+4.0_8*ldos(uu)-3.0_8)*(xa-xb))/jac
	 Miy(4)=(4.0_8*luno(uu)*xa-4.0_8*ldos(uu)*xb)/jac
	 Miy(5)=(4.0_8*ldos(uu)*xb+(4.0_8-4.0_8*luno(uu)-8.0_8*ldos(uu))*xa)/jac
	 Miy(6)=((8.0_8*luno(uu)+4.0_8*ldos(uu)-4.0_8)*xb-4.0_8*luno(uu)*xa)/jac
	 
	 Mp(1)=luno(uu)
	 Mp(2)=ldos(uu)
	 Mp(3)=1.0_8-luno(uu)-ldos(uu)
	 
	  do ui=1,6
  	  do uj=1,3   	 
	  !Posicion 3,1. Caja Bxt.	  
	  Nn(uj,ui)=Nn(uj,ui)+(Mix(ui)*Mp(uj))*ja(uu)*jac 	   
	  !Posicion 3,2. Caja Byt.
	  Nm(uj,ui)=Nm(uj,ui)+(Miy(ui)*Mp(uj))*ja(uu)*jac
	  enddo
	  enddo	   	  
     enddo 
	 call suma(Nn/del,n,2*i,0,3*i,3,6,poi)
	 call suma(Nm/del,n,2*i,i,3*i,3,6,podosi)
 enddo
 close(1)
 deallocate(poi,podosi) 
end

!--------------------------------------------------------------------------------------------------------------------------------------
!Subrutina CAJASDE
!A continuación se calculan coeficientes de las matrices elementales que forman las cajas Dx y Dy.
!'Nn' son las matrices elementales de la caja Dx, 'Nm' son las matrices elementales de la caja Dy. Estas cajas aparecen al discretizar 
!la ecuación de continuidad para las ecuaciones de aguas someras. 
!En 'vt' está el calado solución de la iteración anterior del que dependen estas cajas. Por ello se trata de una caja no lineal del 
!sistema. 'hxy' es el calado evaluado en cualquier punto de integración del elemento, a partir de los valores de calado de la 
!iteración anterior. Aquí sólo se pueden utilizar valores positivos de calado. 
!--------------------------------------------------------------------------------------------------------------------------------------
subroutine cajasde(i,j,x,y,vt,del)
use elemental
integer*4, dimension(:),allocatable::poi,podosi
integer*4 i,j										   
real*8 x(i),y(i),Nn(3,6),Nm(3,6),Mi(6),Mp(3),Mpx(3),Mpy(3),hxy,vt(3*i),del

allocate(poi(3*i),podosi(3*i))
 
39     format(4/,A80)
40     format(6X,6(X,I5))

!Inizialización previa de variables:
!-----------------------------------
luno=(/1.0_8/3.0_8, 0.6_8, 0.2_8, 0.2_8, 0.0_8, 0.0_8, 0.0_8/)		   
ldos=(/1.0_8/3.0_8, 0.2_8, 0.6_8, 0.2_8, 0.0_8, 0.0_8, 0.0_8/)
ja=(/-27.0_8/96.0_8, 25.0_8/96.0_8, 25.0_8/96.0_8, 25.0_8/96.0_8, 0.0_8, 0.0_8, 0.0_8/)
do u=1,3*i
poi(u)=posi(u)
podosi(u)=posdosi(u)
enddo

!Cálculo de los valores para cada elemento
!-----------------------------------------
open(unit=1,file='C:\malla.txt',status='old')
read(1,39)ac
 do u=1,j
     read(1,40) n(1),n(4),n(2),n(5),n(3),n(6)
  
     xa=x(n(1))-x(n(3)) 
	 xb=x(n(2))-x(n(3)) 
	 ya=y(n(1))-y(n(3))
	 yb=y(n(2))-y(n(3))
	 jac=abs(xa*yb-xb*ya)

	 do uu=1,3
     do uuu=1,6
     Nn(uu,uuu)=0.0_8
	 Nm(uu,uuu)=0.0_8
     enddo
     enddo

	 !Integración con cuatro puntos. Polinomios de grado superior a grado 3 quedarían integrados con error.
	 do uu=1,4
	 			   
	 Mi(1)=luno(uu)*(2.0_8*luno(uu)-1.0_8)
	 Mi(2)=ldos(uu)*(2.0_8*ldos(uu)-1.0_8)
	 Mi(3)=(1.0_8-luno(uu)-ldos(uu))*(1.0_8-2.0_8*luno(uu)-2.0_8*ldos(uu))
	 Mi(4)=4.0_8*luno(uu)*ldos(uu)
	 Mi(5)=4.0_8*ldos(uu)*(1.0_8-luno(uu)-ldos(uu))
	 Mi(6)=4.0_8*luno(uu)*(1.0_8-luno(uu)-ldos(uu))

	 Mp(1)=luno(uu)
	 Mp(2)=ldos(uu)
	 Mp(3)=1.0_8-luno(uu)-ldos(uu)

	 Mpx(1)=(y(n(2))-y(n(3)))/jac
	 Mpx(2)=(y(n(3))-y(n(1)))/jac
	 Mpx(3)=(y(n(1))-y(n(2)))/jac
	 
	 Mpy(1)=(x(n(3))-x(n(2)))/jac
	 Mpy(2)=(x(n(1))-x(n(3)))/jac
	 Mpy(3)=(x(n(2))-x(n(1)))/jac 
	 
	 !Se utilizan los calados calculados en la iteracion anterior.
	 hxy=vt(2*i+n(1))*Mp(1)+vt(2*i+n(2))*Mp(2)+vt(2*i+n(3))*Mp(3)
	
	 do ui=1,6
	 do uj=1,3	 
	 !Posicion 3,1. Caja Dx.
	 Nn(uj,ui)=Nn(uj,ui)+(Mpx(uj)*Mi(ui)*hxy)*ja(uu)*jac	 
	 !Posicion 3,2. Caja Dy.
	 Nm(uj,ui)=Nm(uj,ui)+(Mpy(uj)*Mi(ui)*hxy)*ja(uu)*jac
	 enddo
	 enddo		   	  
   enddo
   call suma(-Nn/del,n,2*i,0,3*i,3,6,poi)
   call suma(-Nm/del,n,2*i,i,3*i,3,6,podosi)
 enddo
 close(1)
 deallocate(poi,podosi)
end

!------------------------------------------------------------------------------------------------------------------------------------------
!Subrutina MATRIZNOLINEAL
!A continuación se calculan coeficientes de las matrices elementales que forman la caja C.
!'Nn' son las matrices elementales de esta caja. Estas cajas aparecen al discretizar las ecuaciones dinámicas (ec. para flujo superficial).
!En 'vv' están las velocidades solución de la iteración anterior de las que dependen estas cajas. Por ello se trata de una caja no 
!lineal del sistema.
!'Mu, Mv' son las velocidades en dirección x y en dirección y evaluadas en cualquier punto de integración del elemento, a partir de 
!los valores de velocidad la iteración anterior.
!------------------------------------------------------------------------------------------------------------------------------------------
subroutine matriznolineal(i,j,x,y,vv,del)
use elemental
integer*4, dimension(:),allocatable::poi,podosi
integer*4 i,j
real*8 x(i),y(i),Nn(6,6),Mi(6),Mix(6),Miy(6),Mu,Mv,vv(3*i),del

allocate(poi(3*i),podosi(3*i))

 39     format(4/,A80)
 40     format(6X,6(X,I5))

!Inizialización previa de variables:
!-----------------------------------
luno=(/1.0_8/3.0_8, 0.0597158717_8, 0.4701420641_8, 0.4701420641_8, 0.7974269853_8, 0.1012865073_8, 0.1012865073_8/)		   
ldos=(/1.0_8/3.0_8, 0.4701420641_8, 0.0597158717_8, 0.4701420641_8, 0.1012865073_8, 0.7974269853_8, 0.1012865073_8/)
ja=(/0.1125_8, 0.06619707635_8, 0.06619707635_8, 0.06619707635_8, 0.06296959025_8, 0.06296959025_8, 0.06296959025_8/)  
do u=1,3*i
poi(u)=posi(u)
podosi(u)=posdosi(u)
enddo

!Cálculo de los valores para cada elemento
!-----------------------------------------
open(unit=1,file='C:\malla.txt',status='old')
read(1,39)ac
 do u=1,j
     read(1,40) n(1),n(4),n(2),n(5),n(3),n(6)

	 xa=x(n(1))-x(n(3)) 
	 xb=x(n(2))-x(n(3)) 
	 ya=y(n(1))-y(n(3))
	 yb=y(n(2))-y(n(3))  
	 jac=abs(xa*yb-xb*ya)
	 	 
	 do uu=1,6
     do uuu=1,6
     Nn(uu,uuu)=0.0_8
     enddo
     enddo

	 !Integración gauss de cuarto orden con siete puntos de integración. Polinomios de grado superior a grado 5 quedarían integrados con error.
	 do uu=1,7
	 Mi(1)=luno(uu)*(2.0_8*luno(uu)-1.0_8)
	 Mi(2)=ldos(uu)*(2.0_8*ldos(uu)-1.0_8)
	 Mi(3)=(1.0_8-luno(uu)-ldos(uu))*(1.0_8-2.0_8*luno(uu)-2.0_8*ldos(uu))
	 Mi(4)=4.0_8*luno(uu)*ldos(uu)
	 Mi(5)=4.0_8*ldos(uu)*(1.0_8-luno(uu)-ldos(uu))
	 Mi(6)=4.0_8*luno(uu)*(1.0_8-luno(uu)-ldos(uu))

	 Mix(1)=((4.0_8*luno(uu)-1.0_8)*yb)/jac
	 Mix(2)=((1.0_8-4.0_8*ldos(uu))*ya)/jac
	 Mix(3)=((4.0_8*luno(uu)+4.0_8*ldos(uu)-3.0_8)*(yb-ya))/jac
	 Mix(4)=(4.0_8*ldos(uu)*yb-4.0_8*luno(uu)*ya)/jac
	 Mix(5)=((8.0_8*ldos(uu)+4.0_8*luno(uu)-4.0_8)*ya-4.0_8*ldos(uu)*yb)/jac
	 Mix(6)=((4.0_8-4.0_8*ldos(uu)-8.0_8*luno(uu))*yb+4.0_8*luno(uu)*ya)/jac
	
	 Miy(1)=((1.0_8-4.0_8*luno(uu))*xb)/jac
	 Miy(2)=((4.0_8*ldos(uu)-1.0_8)*xa)/jac
	 Miy(3)=((4.0_8*luno(uu)+4.0_8*ldos(uu)-3.0_8)*(xa-xb))/jac
	 Miy(4)=(4.0_8*luno(uu)*xa-4.0_8*ldos(uu)*xb)/jac
	 Miy(5)=(4.0_8*ldos(uu)*xb+(4.0_8-4.0_8*luno(uu)-8.0_8*ldos(uu))*xa)/jac
	 Miy(6)=((8.0_8*luno(uu)+4.0_8*ldos(uu)-4.0_8)*xb-4.0_8*luno(uu)*xa)/jac

	 !Se calcula el producto escalar
	 Mu=dot_product(vv(n),Mi)   
	 Mv=dot_product(vv(n+i),Mi) 	 	
	 
	  do ui=1,6
  	  do uj=1,6 
	  !Posición 1,1 y 2,2. Caja C.
	  Nn(uj,ui)=Nn(uj,ui)+Mi(uj)*(Mix(ui)*Mu+Miy(ui)*Mv)*ja(uu)*jac
	  enddo
	  enddo	 
	 enddo
	 !La caja C se escribe en dos posiciones diferentes. Además se trata de las posiciones donde se escribió la caja A en la 
	 !subrutina cajasab. De ahí la utilidad de las variables 'poi,podosi'.
	 call suma(Nn/del,n,0,0,3*i,6,6,poi)
	 call suma(Nn/del,n,i,i,3*i,6,6,podosi)
 enddo
 close(1)
 deallocate(poi,podosi)
end
	
!------------------------------------------------------------------------------------------------------------------------------------------
!Subrutina TIMEASU
!A continuación se calculan coeficientes de las matrices elementales que forman la caja M, típicamente llamada matriz de masa.
!'Nn' son las matrices elementales de esta caja. Estas cajas aparecen al discretizar las ecuaciones dinámicas (ec. para flujo superficial).
!Sólo entra aquí en caso de resolver el problema no estacionario. 
!------------------------------------------------------------------------------------------------------------------------------------------
subroutine timeasu(i,j,x,y,At)
use elemental
integer*4, dimension(:),allocatable::poi,podosi
integer*4 i,j				  
real*8 x(i),y(i),Nn(6,6),Mi(6),At

allocate(poi(3*i),podosi(3*i))

 39     format(4/,A80)
 40     format(6X,6(X,I5))

!Inizialización previa de variables:
!-----------------------------------
luno=(/1.0_8/3.0_8, 0.0597158717_8, 0.4701420641_8, 0.4701420641_8, 0.7974269853_8, 0.1012865073_8, 0.1012865073_8/)		   
ldos=(/1.0_8/3.0_8, 0.4701420641_8, 0.0597158717_8, 0.4701420641_8, 0.1012865073_8, 0.7974269853_8, 0.1012865073_8/)
ja=(/0.1125_8, 0.06619707635_8, 0.06619707635_8, 0.06619707635_8, 0.06296959025_8, 0.06296959025_8, 0.06296959025_8/)
do u=1,3*i
poi(u)=posi(u)
podosi(u)=posdosi(u)
enddo

!Cálculo de los valores para cada elemento
!-----------------------------------------
open(unit=1,file='C:\malla.txt',status='old')  
read(1,39)ac	 
 do u=1,j
     read(1,40) n(1),n(4),n(2),n(5),n(3),n(6)
  
	 xa=x(n(1))-x(n(3)) 
	 xb=x(n(2))-x(n(3)) 
	 ya=y(n(1))-y(n(3))
	 yb=y(n(2))-y(n(3))  	 
	 jac=abs(xa*yb-xb*ya)
	 	 
	 do uu=1,6
     do uuu=1,6
     Nn(uu,uuu)=0.0_8
     enddo
     enddo

	 do uu=1,7
	 Mi(1)=luno(uu)*(2.0_8*luno(uu)-1.0_8)
	 Mi(2)=ldos(uu)*(2.0_8*ldos(uu)-1.0_8)
	 Mi(3)=(1.0_8-luno(uu)-ldos(uu))*(1.0_8-2.0_8*luno(uu)-2.0_8*ldos(uu))
	 Mi(4)=4.0_8*luno(uu)*ldos(uu)
	 Mi(5)=4.0_8*ldos(uu)*(1.0_8-luno(uu)-ldos(uu))
	 Mi(6)=4.0_8*luno(uu)*(1.0_8-luno(uu)-ldos(uu))

	  do ui=1,6
  	  do uj=1,6 
	  !Posición 1,1 y 2,2. Caja M.
	  Nn(uj,ui)=Nn(uj,ui)+(Mi(uj)*Mi(ui))*ja(uu)*jac
	  enddo
	  enddo	 
	 enddo
	 call suma(Nn/At,n,0,0,3*i,6,6,poi)
	 call suma(Nn/At,n,i,i,3*i,6,6,podosi)
 enddo
 close(1)
 deallocate(poi,podosi)
end

!------------------------------------------------------------------------------------------------------------------------------
!Subrutina TIMEASD
!A continuación se calculan coeficientes de las matrices elementales que forman la caja N, típicamente llamada matriz de masa.
!'Nm' son las matrices elementales de esta caja. Estas cajas aparecen al discretizar la ecuación de continuidad de las 
!ecuaciones de aguas someras.
!Sólo entra aquí en caso de resolver el problema no estacionario.
!------------------------------------------------------------------------------------------------------------------------------
subroutine timeasd(i,j,x,y,At)
use elemental
integer*4, dimension(:),allocatable::potresi
integer*4 i,j				  
real*8 x(i),y(i),Nm(3,3),Mp(3),At

allocate(potresi(3*i))

 39     format(4/,A80)
 40     format(6X,6(X,I5))

!Inizialización previa de variables:
!-----------------------------------
luno=(/0.5_8, 0.0_8, 0.5_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8/)		   
ldos=(/0.5_8, 0.5_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8/)
ja=(/1.0_8/6.0_8, 1.0_8/6.0_8, 1.0_8/6.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8/)
do u=1,3*i
potresi(u)=postresi(u)
enddo

!Cálculo de los valores para cada elemento
!-----------------------------------------
open(unit=1,file='C:\malla.txt',status='old')  
read(1,39)ac	 
 do u=1,j
     read(1,40) n(1),n(4),n(2),n(5),n(3),n(6)
  
	 xa=x(n(1))-x(n(3)) 
	 xb=x(n(2))-x(n(3)) 
	 ya=y(n(1))-y(n(3))
	 yb=y(n(2))-y(n(3))  	 
	 jac=abs(xa*yb-xb*ya)
	 	 
	 do uu=1,3
     do uuu=1,3
     Nm(uu,uuu)=0.0_8
     enddo
     enddo

	 do uu=1,3
	 Mp(1)=luno(uu)
	 Mp(2)=ldos(uu)
	 Mp(3)=1.0_8-luno(uu)-ldos(uu)

	  do ui=1,3
  	  do uj=1,3
      !Posición 3,3. Caja N.
	  Nm(uj,ui)=Nm(uj,ui)+(Mp(uj)*Mp(ui))*ja(uu)*jac
	  enddo
	  enddo	 
	 enddo
	 call suma(Nm/At,n,2*i,2*i,3*i,3,3,potresi)
 enddo
 close(1)
 deallocate(potresi)
end

!-----------------------------------------------------------------------------------------------------------------------------------------
!Subrutina F
!Con la siguiente subrutina se calcula la fricción por la tensión de fondo a través del número de Manning. 
!Se calculan vectores elementales que irán al término independiente del sistema y que aparecen en las ecuaciones dinámicas de las 
!ecuaciones de aguas someras.
!Con 'yn' se da un rozamiento tal que I=i. Así, independientemente de la velocidad y el calado, se dará un calado normal.
!Esto es equivalente a resolver un sistema con la ecuación de continuidad de aguas someras y las ecuaciones dinámicas de Navier-Stokes 2D.
!-----------------------------------------------------------------------------------------------------------------------------------------
subroutine f(i,j,x,y,z,ma,yn,vv,Nx,Ny)	 
use elemental
integer*4 i,j	 
real*8 x(i),y(i),Nx(i),Ny(i),s,maning,ma(i),vv(3*i),Mi(6),Mu,Mv,Mui,Mvi,Mhi,dist,a,b,c
real*8 z(i),Mix(6),Miy(6),Mpx,Mpy  
character yn*2
			 
39     format(4/,A80)
40     format(6X,6(X,I5))
		  													 
!Inizialización previa de variables:
!-----------------------------------
luno=(/1.0_8/3.0_8, 0.0597158717_8, 0.4701420641_8, 0.4701420641_8, 0.7974269853_8, 0.1012865073_8, 0.1012865073_8/)		   
ldos=(/1.0_8/3.0_8, 0.4701420641_8, 0.0597158717_8, 0.4701420641_8, 0.1012865073_8, 0.7974269853_8, 0.1012865073_8/)
ja=(/0.1125_8, 0.06619707635_8, 0.06619707635_8, 0.06619707635_8, 0.06296959025_8, 0.06296959025_8, 0.06296959025_8/)    		   

!Cálculo de los valores para cada elemento
!-----------------------------------------
open(unit=1,file='C:\malla.txt',status='old')
read(1,39)ac
 do u=1,j
 read(1,40) n(1),n(4),n(2),n(5),n(3),n(6)

	 xa=x(n(1))-x(n(3)) 
	 xb=x(n(2))-x(n(3)) 
	 ya=y(n(1))-y(n(3))
	 yb=y(n(2))-y(n(3))  
	 jac=abs(xa*yb-xb*ya)	

	 !Cálculo del módulo de la velocidad en cada nodo esquina.
	 a=sqrt(vv(n(1))**2+vv(i+n(1))**2)
	 b=sqrt(vv(n(2))**2+vv(i+n(2))**2)
	 c=sqrt(vv(n(3))**2+vv(i+n(3))**2)
	 !La distancia elemental se estima como el diámetro del círculo de área igual a la del elemento.
	 dist=sqrt((2.0_8*jac)/(3.14159_8))
	 
	 !Valores medios de velocidad, calado y manning en el elemento (valores en el centro del elemento).
	 Mui=-(vv(n(1))+vv(n(2))+vv(n(3)))/9.0_8+(vv(n(4))+vv(n(5))+vv(n(6)))*4.0_8/9.0_8   
	 Mvi=-(vv(n(1)+i)+vv(n(2)+i)+vv(n(3)+i))/9.0_8+(vv(n(4)+i)+vv(n(5)+i)+vv(n(6)+i))*4.0_8/9.0_8
	 Mhi=(vv(n(1)+2*i)+vv(n(2)+2*i)+vv(n(3)+2*i))/3.0_8
	 maning=(ma(n(1))+ma(n(2))+ma(n(3)))/3.0_8
	 
	 !Se da un coeficiente de Manning atendiendo al tipo de elemento.
	 if (((a.eq.0.0).and.(b.eq.0.0)).or.((a.eq.0.0).and.(c.eq.0.0)).or.((b.eq.0.0).and.(c.eq.0.0))) then
	  !Elementos pegados a un contorno con velocidades nulas. Se consideran sólo elementos con dos nodos apoyados en la pared. 
	  !El perímetro mojado tiene en cuenta la pared lateral y el radio hidráulico es función de la distancia. 
	  s=((maning**2.0_8))/(((dist*Mhi)/(dist+Mhi))**(4.0_8/3.0_8))	  
	 else 
	  !Se consideran elementos con un nodo apoyado en la pared o con ningún nodo apoyado.
	  !En aquéllos con un nodo apoyado se tiene la pared en un punto infinitesimal del elemento. Ésta ya se considera en elementos con dos 
	  !nodos apoyados. El radio hidráulico se estima como el calado medio.
	  s=((maning**2.0_8))/(Mhi**(4.0_8/3.0_8)) 
	 endif
	  
	 do uu=1,7	      
	 Mi(1)=luno(uu)*(2.0_8*luno(uu)-1.0_8)
	 Mi(2)=ldos(uu)*(2.0_8*ldos(uu)-1.0_8)
	 Mi(3)=(1.0_8-luno(uu)-ldos(uu))*(1.0_8-2.0_8*luno(uu)-2.0_8*ldos(uu))
	 Mi(4)=4.0_8*luno(uu)*ldos(uu)
	 Mi(5)=4.0_8*ldos(uu)*(1.0_8-luno(uu)-ldos(uu))
	 Mi(6)=4.0_8*luno(uu)*(1.0_8-luno(uu)-ldos(uu))

	 !Procedimiento con poco sentido físico.
	 !Se obliga a que Ix=ix, Iy=iy (con i cambiado de signo, ix será positivo si z baja a medida que crecen los valores de x).
	 !Se ha comprobado que efectivamente con esto la altura de la lámina de agua se ajusta a la cota del terreno (calado normal).
	 if (yn.eq.'si') then
	 Mix(1)=((4.0_8*luno(uu)-1.0_8)*yb)/jac
	 Mix(2)=((1.0_8-4.0_8*ldos(uu))*ya)/jac
	 Mix(3)=((4.0_8*luno(uu)+4.0_8*ldos(uu)-3.0_8)*(yb-ya))/jac
	 Mix(4)=(4.0_8*ldos(uu)*yb-4.0_8*luno(uu)*ya)/jac
	 Mix(5)=((8.0_8*ldos(uu)+4.0_8*luno(uu)-4.0_8)*ya-4.0_8*ldos(uu)*yb)/jac
	 Mix(6)=((4.0_8-4.0_8*ldos(uu)-8.0_8*luno(uu))*yb+4.0_8*luno(uu)*ya)/jac
	
	 Miy(1)=((1.0_8-4.0_8*luno(uu))*xb)/jac
	 Miy(2)=((4.0_8*ldos(uu)-1.0_8)*xa)/jac
	 Miy(3)=((4.0_8*luno(uu)+4.0_8*ldos(uu)-3.0_8)*(xa-xb))/jac
	 Miy(4)=(4.0_8*luno(uu)*xa-4.0_8*ldos(uu)*xb)/jac
	 Miy(5)=(4.0_8*ldos(uu)*xb+(4.0_8-4.0_8*luno(uu)-8.0_8*ldos(uu))*xa)/jac
	 Miy(6)=((8.0_8*luno(uu)+4.0_8*ldos(uu)-4.0_8)*xb-4.0_8*luno(uu)*xa)/jac

	  Mpx=dot_product(z(n),Mix)	  
	  Mpy=dot_product(z(n),Miy)
	  do uuu=1,6
	  Nx(n(uuu))=Nx(n(uuu))-Mi(uuu)*Mpx*jac*ja(uu) 	 
	  Ny(n(uuu))=Ny(n(uuu))-Mi(uuu)*Mpy*jac*ja(uu) 
	  enddo
	 
	 !Procedimiento para calcular con las ecuaciones de aguas someras. 
	 !Se calcula la pendiente de fricción a partir de las variables de velocidad y calado.
	 else
	  Mu=dot_product(vv(n),Mi)   
	  Mv=dot_product(vv(n+i),Mi) 	  
	  do uuu=1,6
	  Nx(n(uuu))=Nx(n(uuu))+(Mi(uuu)*Mu*sqrt(Mui**2+Mvi**2)*s)*ja(uu)*jac 	 
	  Ny(n(uuu))=Ny(n(uuu))+(Mi(uuu)*Mv*sqrt(Mui**2+Mvi**2)*s)*ja(uu)*jac 
	  enddo	
	 endif
	 enddo
 enddo 
 close(1)
 end

!------------------------------------------------------------------------------------------------------------------------
!Subrutina JACOB
!Esta subrutina calcula los términos que faltan de la matriz jacobiana del sistema para las ecuaciones de aguas someras.
!También se dan condiciones de contorno sobre el sistema como en la subrutina reducciondelsistema.
!Así, se tendrá preparado el sistema a resolver para el método de Newton.
!------------------------------------------------------------------------------------------------------------------------
subroutine jacob (i,j,x,y,vv,nu,man,vector,v,vic,vb,cc,ten,bcg,ndim)
use elemental
use allocatacion
integer*4, dimension(:),allocatable::poi,podosi,potresi,ck
integer*4 i,j,ndim,v(i),cc,inc,w,k
real*8 x(i),y(i),Ncua(6,6),Ncub(6,6),Ncva(6,6),Ncvb(6,6),Nbhx(6,3),Nbhy(6,3),Nbuh(3,6),Nbvh(3,6),Nbu(3,3),Nbv(3,3)
real*8 Ndyx(3,3),Nfx(6,6),Nfy(6,6),Nbnux(6,6),Mi(6),Mix(6),Miy(6),Mp(3),Mpx(3),Mpy(3),vv(3*i),vector(3*i),vic(3*i)
real*8 Mu,Mv,Mux,Mvx,Muy,Mvy,Mui,Mvi,Mhi,maning,man(i),amn,bmn,akm,bkm,ank,bnk,a,b,c,dist,s,ma,mb,mc,md,me,nu,vb(3*i)
character ten*2,bcg*2

allocate(poi(3*i),podosi(3*i),potresi(3*i),ck(3*i))

 39     format(4/,A80)										  
 40     format(6X,6(X,I5))

!Inizialización previa de variables:
!-----------------------------------
!Se utilizarán dos integraciones diferentes en el dominio, una con 'luno,ldos,ja' y otra con 'lunot,ldost,jat'.
luno=(/1.0_8/3.0_8, 0.6_8, 0.2_8, 0.2_8, 0.0_8, 0.0_8, 0.0_8/)		   
ldos=(/1.0_8/3.0_8, 0.2_8, 0.6_8, 0.2_8, 0.0_8, 0.0_8, 0.0_8/)
ja=(/-27.0_8/96.0_8, 25.0_8/96.0_8, 25.0_8/96.0_8, 25.0_8/96.0_8, 0.0_8, 0.0_8, 0.0_8/)
lunot=(/1.0_8/3.0_8, 0.0597158717_8, 0.4701420641_8, 0.4701420641_8, 0.7974269853_8, 0.1012865073_8, 0.1012865073_8/)		   
ldost=(/1.0_8/3.0_8, 0.4701420641_8, 0.0597158717_8, 0.4701420641_8, 0.1012865073_8, 0.7974269853_8, 0.1012865073_8/)
jat=(/0.1125_8, 0.06619707635_8, 0.06619707635_8, 0.06619707635_8, 0.06296959025_8, 0.06296959025_8, 0.06296959025_8/)
do u=1,3*i
poi(u)=posi(u)
podosi(u)=posdosi(u)
potresi(u)=postresi(u)
enddo
ma=9.0_8
mb=12.0_8
mc=8.0_8
md=3.0_8
me=4.0_8

open(unit=1,file='C:\malla.txt',status='old')
read(1,39)ac
 do u=1,j
     read(1,40) n(1),n(4),n(2),n(5),n(3),n(6)

	 xa=x(n(1))-x(n(3)) 
	 xb=x(n(2))-x(n(3)) 
	 ya=y(n(1))-y(n(3))
	 yb=y(n(2))-y(n(3))  
	 jac=abs(xa*yb-xb*ya)
	 
	 !Cálculo de parámetros para el cálculo del término correspondiente a la pendiente de fricción, de un modo similar al de la subrutina f.
	 !Se considerará sólo el caso yn='no' indicado en la subrutina f.
	 a=sqrt(vv(n(1))**2+vv(i+n(1))**2)
	 b=sqrt(vv(n(2))**2+vv(i+n(2))**2)
	 c=sqrt(vv(n(3))**2+vv(i+n(3))**2)
	 dist=sqrt((2.0_8*jac)/(3.14159_8))
	 Mui=-(vv(n(1))+vv(n(2))+vv(n(3)))/9.0_8+(vv(n(4))+vv(n(5))+vv(n(6)))*4.0_8/9.0_8   
	 Mvi=-(vv(n(1)+i)+vv(n(2)+i)+vv(n(3)+i))/9.0_8+(vv(n(4)+i)+vv(n(5)+i)+vv(n(6)+i))*4.0_8/9.0_8
	 Mhi=(vv(n(1)+2*i)+vv(n(2)+2*i)+vv(n(3)+2*i))/3.0_8
	 maning=(man(n(1))+man(n(2))+man(n(3)))/3.0_8
	  if (((a.eq.0.0).and.(b.eq.0.0)).or.((a.eq.0.0).and.(c.eq.0.0)).or.((b.eq.0.0).and.(c.eq.0.0))) then
	   s=((maning**2.0_8))/(((dist*Mhi)/(dist+Mhi))**(4.0_8/3.0_8))	  
	  else 
	   s=((maning**2.0_8))/(Mhi**(4.0_8/3.0_8)) 
	  endif

	 do uu=1,6
     do uuu=1,6
	 Ncua(uu,uuu)=0.0_8
	 Ncub(uu,uuu)=0.0_8
	 Ncva(uu,uuu)=0.0_8
	 Ncvb(uu,uuu)=0.0_8 
	 Nfx(uu,uuu)=0.0_8
	 Nfy(uu,uuu)=0.0_8	
     enddo
     enddo
	 do uu=1,3
     do uuu=1,3
	 Ndyx(uu,uuu)=0.0_8	
     enddo
     enddo

	 call direccion(x,y,i,n(2),n(3),n(1),amn,bmn)
	 call direccion(x,y,i,n(1),n(2),n(3),akm,bkm)
	 call direccion(x,y,i,n(3),n(1),n(2),ank,bnk)
	 
	 !Cajas formadas a partir de integrales de contorno
	 !-------------------------------------------------
	 !Derivadas de las integrales de contorno de los términos de presión para el método de Newton.
	 Nbhx(1,1)=(akm*abs(ya-yb)+ank*abs(ya))/6.0_8	 
	 Nbhx(2,2)=(akm*abs(ya-yb)+amn*abs(yb))/6.0_8	 
	 Nbhx(3,3)=(amn*abs(yb)+ank*abs(ya))/6.0_8
	 Nbhx(4,1)=akm*abs(ya-yb)/3.0_8
	 Nbhx(4,2)=akm*abs(ya-yb)/3.0_8
	 Nbhx(5,2)=amn*abs(yb)/3.0_8
	 Nbhx(5,3)=amn*abs(yb)/3.0_8
	 Nbhx(6,1)=ank*abs(ya)/3.0_8
	 Nbhx(6,3)=ank*abs(ya)/3.0_8
	 
	 Nbhy(1,1)=(bkm*abs(xa-xb)+bnk*abs(xa))/6.0_8	 
	 Nbhy(2,2)=(bkm*abs(xa-xb)+bmn*abs(xb))/6.0_8	 
	 Nbhy(3,3)=(bmn*abs(xb)+bnk*abs(xa))/6.0_8
	 Nbhy(4,1)=bkm*abs(xa-xb)/3.0_8
	 Nbhy(4,2)=bkm*abs(xa-xb)/3.0_8
	 Nbhy(5,2)=bmn*abs(xb)/3.0_8
	 Nbhy(5,3)=bmn*abs(xb)/3.0_8
	 Nbhy(6,1)=bnk*abs(xa)/3.0_8
	 Nbhy(6,3)=bnk*abs(xa)/3.0_8

	 !Derivadas de las integrales de contorno de los términos viscosos para el método de Newton. 
	 if (ten.eq.'si') then 
	 !Se consideran si también se consideran estas integrales de contorno
	 Nbnux(1,1)=md*(yb*ank*abs(ya)+ya*akm*abs(ya-yb)-xb*bnk*abs(xa)-xa*bkm*abs(xa-xb))/jac
	 Nbnux(1,2)=(ya*akm*abs(ya-yb)-xa*bkm*abs(xa-xb))/jac
	 Nbnux(1,3)=(yb*ank*abs(ya)-xb*bnk*abs(xa))/jac
     Nbnux(1,4)=me*(-ya*akm*abs(ya-yb)+xa*bkm*abs(xa-xb))/jac
     Nbnux(1,6)=me*(-yb*ank*abs(ya)+xb*bnk*abs(xa))/jac
	 Nbnux(2,1)=(-ya*akm*abs(ya-yb)+xa*bkm*abs(xa-xb))/jac
	 Nbnux(2,2)=md*(-ya*akm*abs(ya-yb)-(ya-yb)*amn*abs(yb)+xa*bkm*abs(xa-xb)+(xa-xb)*bmn*abs(xb))/jac
	 Nbnux(2,3)=(-(ya-yb)*amn*abs(yb)+(xa-xb)*bmn*abs(xb))/jac
	 Nbnux(2,4)=me*(ya*akm*abs(ya-yb)-xa*bkm*abs(xa-xb))/jac
	 Nbnux(2,5)=me*((ya-yb)*amn*abs(yb)-(xa-xb)*bmn*abs(xb))/jac
	 Nbnux(3,1)=(-yb*ank*abs(ya)+xb*bnk*abs(xa))/jac
	 Nbnux(3,2)=((ya-yb)*amn*abs(yb)-(xa-xb)*bmn*abs(xb))/jac
	 Nbnux(3,3)=md*(-yb*ank*abs(ya)+(ya-yb)*amn*abs(yb)+xb*bnk*abs(xa)-(xa-xb)*bmn*abs(xb))/jac
	 Nbnux(3,5)=me*(-(ya-yb)*amn*abs(yb)+(xa-xb)*bmn*abs(xb))/jac
	 Nbnux(3,6)=me*(yb*ank*abs(ya)-xb*bnk*abs(xa))/jac
     Nbnux(4,1)=me*(ya*akm*abs(ya-yb)-xa*bkm*abs(xa-xb))/jac
     Nbnux(4,2)=me*(-ya*akm*abs(ya-yb)+xa*bkm*abs(xa-xb))/jac
	 Nbnux(5,2)=me*(-(ya-yb)*amn*abs(yb)+(xa-xb)*bmn*abs(xb))/jac
	 Nbnux(5,3)=me*((ya-yb)*amn*abs(yb)-(xa-xb)*bmn*abs(xb))/jac
	 Nbnux(6,1)=me*(yb*ank*abs(ya)-xb*bnk*abs(xa))/jac
	 Nbnux(6,3)=me*(-yb*ank*abs(ya)+xb*bnk*abs(xa))/jac
	 else
	 do uu=1,6
     do uuu=1,6
	 Nbnux(uu,uuu)=0.0_8	
     enddo
     enddo
	 endif

	 !Derivadas de las integrales de contorno de la ecuación de continuidad para el método de Newton.
	 Nbuh(1,1)=(ma*vv(2*i+n(1))+vv(2*i+n(2)))*akm*abs(ya-yb)+(ma*vv(2*i+n(1))+vv(2*i+n(3)))*ank*abs(ya)
	 Nbuh(1,2)=(-vv(2*i+n(1))+vv(2*i+n(2)))*akm*abs(ya-yb)
	 Nbuh(1,3)=(-vv(2*i+n(1))+vv(2*i+n(3)))*ank*abs(ya)
	 Nbuh(1,4)=(mb*vv(2*i+n(1))+mc*vv(2*i+n(2)))*akm*abs(ya-yb)
	 Nbuh(1,6)=(mb*vv(2*i+n(1))+mc*vv(2*i+n(3)))*ank*abs(ya)
	 Nbuh(2,1)=(vv(2*i+n(1))-vv(2*i+n(2)))*akm*abs(ya-yb)
	 Nbuh(2,2)=(vv(2*i+n(1))+ma*vv(2*i+n(2)))*akm*abs(ya-yb)+(ma*vv(2*i+n(2))+vv(2*i+n(3)))*amn*abs(yb)
	 Nbuh(2,3)=(-vv(2*i+n(2))+vv(2*i+n(3)))*amn*abs(yb)
	 Nbuh(2,4)=(mc*vv(2*i+n(1))+mb*vv(2*i+n(2)))*akm*abs(ya-yb)
	 Nbuh(2,5)=(mb*vv(2*i+n(2))+mc*vv(2*i+n(3)))*amn*abs(yb)	 
	 Nbuh(3,1)=(vv(2*i+n(1))-vv(2*i+n(3)))*ank*abs(ya)
	 Nbuh(3,2)=(vv(2*i+n(2))-vv(2*i+n(3)))*amn*abs(yb)
	 Nbuh(3,3)=(vv(2*i+n(2))+ma*vv(2*i+n(3)))*amn*abs(yb)+(vv(2*i+n(1))+ma*vv(2*i+n(3)))*ank*abs(ya)
	 Nbuh(3,5)=(mc*vv(2*i+n(2))+mb*vv(2*i+n(3)))*amn*abs(yb)
	 Nbuh(3,6)=(mc*vv(2*i+n(1))+mb*vv(2*i+n(3)))*ank*abs(ya)
	 
	 Nbvh(1,1)=(ma*vv(2*i+n(1))+vv(2*i+n(2)))*bkm*abs(xa-xb)+(ma*vv(2*i+n(1))+vv(2*i+n(3)))*bnk*abs(xa)
	 Nbvh(1,2)=(-vv(2*i+n(1))+vv(2*i+n(2)))*bkm*abs(xa-xb)
	 Nbvh(1,3)=(-vv(2*i+n(1))+vv(2*i+n(3)))*bnk*abs(xa)
	 Nbvh(1,4)=(mb*vv(2*i+n(1))+mc*vv(2*i+n(2)))*bkm*abs(xa-xb)
	 Nbvh(1,6)=(mb*vv(2*i+n(1))+mc*vv(2*i+n(3)))*bnk*abs(xa)
	 Nbvh(2,1)=(vv(2*i+n(1))-vv(2*i+n(2)))*bkm*abs(xa-xb)
	 Nbvh(2,2)=(vv(2*i+n(1))+ma*vv(2*i+n(2)))*bkm*abs(xa-xb)+(ma*vv(2*i+n(2))+vv(2*i+n(3)))*bmn*abs(xb)
	 Nbvh(2,3)=(-vv(2*i+n(2))+vv(2*i+n(3)))*bmn*abs(xb)
	 Nbvh(2,4)=(mc*vv(2*i+n(1))+mb*vv(2*i+n(2)))*bkm*abs(xa-xb)
	 Nbvh(2,5)=(mb*vv(2*i+n(2))+mc*vv(2*i+n(3)))*bmn*abs(xb)	 
	 Nbvh(3,1)=(vv(2*i+n(1))-vv(2*i+n(3)))*bnk*abs(xa)
	 Nbvh(3,2)=(vv(2*i+n(2))-vv(2*i+n(3)))*bmn*abs(xb)
	 Nbvh(3,3)=(vv(2*i+n(2))+ma*vv(2*i+n(3)))*bmn*abs(xb)+(vv(2*i+n(1))+ma*vv(2*i+n(3)))*bnk*abs(xa)
	 Nbvh(3,5)=(mc*vv(2*i+n(2))+mb*vv(2*i+n(3)))*bmn*abs(xb)
	 Nbvh(3,6)=(mc*vv(2*i+n(1))+mb*vv(2*i+n(3)))*bnk*abs(xa)

	 Nbu(1,1)=(ma*vv(n(1))-vv(n(2))+mb*vv(n(4)))*akm*abs(ya-yb)+(ma*vv(n(1))-vv(n(3))+mb*vv(n(6)))*ank*abs(ya)
	 Nbu(1,2)=(vv(n(1))+vv(n(2))+mc*vv(n(4)))*akm*abs(ya-yb)
	 Nbu(1,3)=(vv(n(1))+vv(n(3))+mc*vv(n(6)))*ank*abs(ya)
	 Nbu(2,1)=(vv(n(1))+vv(n(2))+mc*vv(n(4)))*akm*abs(ya-yb)
	 Nbu(2,2)=(-vv(n(1))+ma*vv(n(2))+mb*vv(n(4)))*akm*abs(ya-yb)+(ma*vv(n(2))-vv(n(3))+mb*vv(n(5)))*amn*abs(yb)
	 Nbu(2,3)=(vv(n(2))+vv(n(3))+mc*vv(n(5)))*amn*abs(yb)
	 Nbu(3,1)=(vv(n(1))+vv(n(3))+mc*vv(n(6)))*ank*abs(ya)
	 Nbu(3,2)=(vv(n(2))+vv(n(3))+mc*vv(n(5)))*amn*abs(yb)
	 Nbu(3,3)=(-vv(n(2))+ma*vv(n(3))+mb*vv(n(5)))*amn*abs(yb)+(-vv(n(1))+ma*vv(n(3))+mb*vv(n(6)))*ank*abs(ya)

	 Nbv(1,1)=(ma*vv(i+n(1))-vv(i+n(2))+mb*vv(i+n(4)))*bkm*abs(xa-xb)+(ma*vv(i+n(1))-vv(i+n(3))+mb*vv(i+n(6)))*bnk*abs(xa)
	 Nbv(1,2)=(vv(i+n(1))+vv(i+n(2))+mc*vv(i+n(4)))*bkm*abs(xa-xb)
	 Nbv(1,3)=(vv(i+n(1))+vv(i+n(3))+mc*vv(i+n(6)))*bnk*abs(xa)
	 Nbv(2,1)=(vv(i+n(1))+vv(i+n(2))+mc*vv(i+n(4)))*bkm*abs(xa-xb)
	 Nbv(2,2)=(-vv(i+n(1))+ma*vv(i+n(2))+mb*vv(i+n(4)))*bkm*abs(xa-xb)+(ma*vv(i+n(2))-vv(i+n(3))+mb*vv(i+n(5)))*bmn*abs(xb)
	 Nbv(2,3)=(vv(i+n(2))+vv(i+n(3))+mc*vv(i+n(5)))*bmn*abs(xb)
	 Nbv(3,1)=(vv(i+n(1))+vv(i+n(3))+mc*vv(i+n(6)))*bnk*abs(xa)
	 Nbv(3,2)=(vv(i+n(2))+vv(i+n(3))+mc*vv(i+n(5)))*bmn*abs(xb)
	 Nbv(3,3)=(-vv(i+n(2))+ma*vv(i+n(3))+mb*vv(i+n(5)))*bmn*abs(xb)+(-vv(i+n(1))+ma*vv(i+n(3))+mb*vv(i+n(6)))*bnk*abs(xa)

    !Cálculo de los valores de otras matrices elementales (4 puntos de integración)
    !------------------------------------------------------------------------------
	do uu=1,4			   
	 Mi(1)=luno(uu)*(2.0_8*luno(uu)-1.0_8)
	 Mi(2)=ldos(uu)*(2.0_8*ldos(uu)-1.0_8)
	 Mi(3)=(1.0_8-luno(uu)-ldos(uu))*(1.0_8-2.0_8*luno(uu)-2.0_8*ldos(uu))
	 Mi(4)=4.0_8*luno(uu)*ldos(uu)
	 Mi(5)=4.0_8*ldos(uu)*(1.0_8-luno(uu)-ldos(uu))
	 Mi(6)=4.0_8*luno(uu)*(1.0_8-luno(uu)-ldos(uu))

	 Mp(1)=luno(uu)
	 Mp(2)=ldos(uu)
	 Mp(3)=1.0_8-luno(uu)-ldos(uu)

	 Mpx(1)=(y(n(2))-y(n(3)))/jac
	 Mpx(2)=(y(n(3))-y(n(1)))/jac
	 Mpx(3)=(y(n(1))-y(n(2)))/jac
	 
	 Mpy(1)=(x(n(3))-x(n(2)))/jac
	 Mpy(2)=(x(n(1))-x(n(3)))/jac
	 Mpy(3)=(x(n(2))-x(n(1)))/jac

	 Mu=dot_product(vv(n),Mi)   
	 Mv=dot_product(vv(n+i),Mi)

	 do ui=1,3					   
  	 do uj=1,3
     !Posición 3,3. Caja Dyx.
	 Ndyx(uj,ui)=Ndyx(uj,ui)+(Mpx(uj)*Mp(ui)*Mu+Mpy(uj)*Mp(ui)*Mv)*ja(uu)*jac 
	 enddo
	 enddo 		   	  
    enddo

	!Cálculo de los valores de otras matrices elementales (7 puntos de integración)
    !------------------------------------------------------------------------------
    do uu=1,7
	 Mi(1)=lunot(uu)*(2.0_8*lunot(uu)-1.0_8)
	 Mi(2)=ldost(uu)*(2.0_8*ldost(uu)-1.0_8)
	 Mi(3)=(1.0_8-lunot(uu)-ldost(uu))*(1.0_8-2.0_8*lunot(uu)-2.0_8*ldost(uu))
	 Mi(4)=4.0_8*lunot(uu)*ldost(uu)
	 Mi(5)=4.0_8*ldost(uu)*(1.0_8-lunot(uu)-ldost(uu))
	 Mi(6)=4.0_8*lunot(uu)*(1.0_8-lunot(uu)-ldost(uu))

	 Mix(1)=((4.0_8*lunot(uu)-1.0_8)*yb)/jac
	 Mix(2)=((1.0_8-4.0_8*ldost(uu))*ya)/jac
	 Mix(3)=((4.0_8*lunot(uu)+4.0_8*ldost(uu)-3.0_8)*(yb-ya))/jac
	 Mix(4)=(4.0_8*ldost(uu)*yb-4.0_8*lunot(uu)*ya)/jac
	 Mix(5)=((8.0_8*ldost(uu)+4.0_8*lunot(uu)-4.0_8)*ya-4.0_8*ldost(uu)*yb)/jac
	 Mix(6)=((4.0_8-4.0_8*ldost(uu)-8.0_8*lunot(uu))*yb+4.0_8*lunot(uu)*ya)/jac
	
	 Miy(1)=((1.0_8-4.0_8*lunot(uu))*xb)/jac
	 Miy(2)=((4.0_8*ldost(uu)-1.0_8)*xa)/jac
	 Miy(3)=((4.0_8*lunot(uu)+4.0_8*ldost(uu)-3.0_8)*(xa-xb))/jac
	 Miy(4)=(4.0_8*lunot(uu)*xa-4.0_8*ldost(uu)*xb)/jac
	 Miy(5)=(4.0_8*ldost(uu)*xb+(4.0_8-4.0_8*lunot(uu)-8.0_8*ldost(uu))*xa)/jac
	 Miy(6)=((8.0_8*lunot(uu)+4.0_8*ldost(uu)-4.0_8)*xb-4.0_8*lunot(uu)*xa)/jac

	 Mux=dot_product(vv(n),Mix)   
	 Mvx=dot_product(vv(n+i),Mix)
	 Muy=dot_product(vv(n),Miy)   
	 Mvy=dot_product(vv(n+i),Miy)
	 
     do ui=1,6
  	 do uj=1,6 	         
	 !Posición 1,1. Cajas DCUA y término de fricción (Manning).
	 Ncua(uj,ui)=Ncua(uj,ui)+Mi(uj)*(Mi(ui)*Mux)*jat(uu)*jac					 
	 Nfx(uj,ui)=Nfx(uj,ui)+(Mi(uj)*Mi(ui)*sqrt(Mui**2+Mvi**2)*s)*jat(uu)*jac  		         
	 !Posición 1,2. Cajas DCVA. 
	 Ncva(uj,ui)=Ncva(uj,ui)+Mi(uj)*(Mi(ui)*Muy)*jat(uu)*jac			         
	 !Posición 2,1. Cajas DCUB. 
	 Ncub(uj,ui)=Ncub(uj,ui)+Mi(uj)*(Mi(ui)*Mvx)*jat(uu)*jac			         
	 !Posición 2,2. Cajas DCVB y término de fricción (Manning).
	 Ncvb(uj,ui)=Ncvb(uj,ui)+Mi(uj)*(Mi(ui)*Mvy)*jat(uu)*jac					 
	 Nfy(uj,ui)=Nfy(uj,ui)+(Mi(uj)*Mi(ui)*sqrt(Mui**2+Mvi**2)*s)*jat(uu)*jac	 	         
	 enddo
	 enddo	   	  
    enddo
	 call suma(Ncua+Nfx*9.81_8-Nbnux*nu/6.0_8,n,0,0,3*i,6,6,poi)
	 call suma(Ncva,n,0,i,3*i,6,6,podosi)	 
	 call suma(Nbhx*9.81_8,n,0,2*i,3*i,6,3,potresi)	 
	 call suma(Ncub,n,i,0,3*i,6,6,poi)
	 call suma(Ncvb+Nfy*9.81_8-Nbnux*nu/6.0_8,n,i,i,3*i,6,6,podosi)	 
	 call suma(Nbhy*9.81_8,n,i,2*i,3*i,6,3,potresi)
	 call suma(Nbuh/60.0_8,n,2*i,0,3*i,3,6,poi)
     call suma(Nbvh/60.0_8,n,2*i,i,3*i,3,6,podosi)	 
	 call suma(-Ndyx+Nbu/60.0_8+Nbv/60.0_8,n,2*i,2*i,3*i,3,3,potresi)	 
 enddo
 close(1)
!A continuación se siguen los procesos que se siguen en la subrutina reducciondelsistema.

!Escritura de la matriz del sistema en formato MSR (Modified Sparse Row)
!-----------------------------------------------------------------------
!Se trata de un almacenamiento indexado por filas con dos vectores.
inc=3*i
!Generación de vectores con formato MSR.
call orden (inc,k)
ndim=ndim-k

!Imposición de condiciones de contorno sobre el sistema
!------------------------------------------------------ 
cc=0
!Sólo se almacenan en 'vic' CC nulas en las filas y columnas a eliminar por tenerse una discretización menor.
do u=1,i
 if (v(u).eq.1) then
 vic(u+inc-i)=0.0_8
 endif
enddo 

!Trozo diferente al de la subrutina reduccióndelsistema. Aquí no es necesario modificar el término independiente con las CC almacenadas en 'vic' 
!como en la subrutina reduccióndelsistema, pues se introducirán valores nulos como valores conocidos de las incógnitas. 
!Por tanto, no es necesario calcular vi y simplemente se reduce el orden del término independiente (se eliminan filas).   
do u=1,inc
 !El sistema se reduce en aquellos nodos donde haya CC (valor distinto de raíz de dos) no nula o CC nula.
 if (vic(u).eq.sqrt(2.0_8)) then 
 cc=cc+1
 vector(cc)=vector(u) 
 vb(cc)=vb(u)
 endif
enddo

!Se reduce el orden de la matriz del sistema (se eliminan filas y columnas). Se mantiene el formato MSR durante este proceso.
call reduccion(ndim,vic,inc)	  

write(6,*) ' ' 
write(6,*)'Longitud del vector matriz dispersa:',ndim

!Escritura de la matriz del sistema en formato CSC (Compressed Sparse Column)
!----------------------------------------------------------------------------
!Se trata de un almacenamiento indexado por columnas con tres vectores.
if (bcg.eq.'no') then
 !Se deja el formato MSR si se va a calcular con precondicionador diagonal (requerimientos de la programación del solver).
 !Vectores 'ita,sa' con coeficientes hasta 'ndim' y con formato MSR. Aunque su dimensión es mucho mayor se utilizarán estos vectores y sólo 
 !esos coeficientes evitando generar nuevos vectores de dimensión 'ndim' (copiando 'ita,sa' a 'ila,la' tras allocate(ila(ndim),la(ndim))).  
 !Los vectores 'cia,ca' no son necesarios.
 deallocate(cia,ca)
else
 !Se pasa a formato CSC si se va a calcular con precondicionador LU (requerimientos de la programación del solver).
 allocate(cja(cc+1))
 !Inicialización de variables.
 do u=1,3*i
 ck(u)=0
 enddo
 cja(1)=1
 do u=1,cc
 cja(u+1)=0
 enddo
 !Cálculo del número de coeficientes por columna.
 do u=cc+2,ndim
 cja(ita(u)+1)=cja(ita(u)+1)+1
 enddo
 !Configuración final de 'cja' (considerando el coeficiente de la diagonal) y escritura de la diagonal en 'cia,ca'.
 !Además se utiliza el contador 'ck'.
 do u=1,cc						
 cja(u+1)=cja(u+1)+cja(u)+1
 ca(cja(u))=sa(u)
 cia(cja(u))=u
 ck(u)=cja(u)+1
 enddo
 !Escritura del resto de coeficientes en 'cia,ca'.
 do u=1,cc										   
  do w=ita(u),ita(u+1)-1
  !En 'w' está la columna y en 'ck(w)' la posición preparada para ese coeficiente.
  cia(ck(ita(w)))=u
  ca(ck(ita(w)))=sa(w)
  ck(ita(w))=ck(ita(w))+1
  enddo
 enddo
 !Vectores 'cia,ca' con coeficientes hasta 'ndim'-1 y con formato CSC. Aunque su dimensión es mucho mayor se utilizarán estos vectores y sólo esos 
 !coeficientes evitando generar nuevos vectores de dimensión 'ndim'-1 (copiando 'ita,sa' a 'dia,da' tras allocate(dia(ndim-1),da(ndim-1),cja(cc+1))). 
 !Los vectores 'ita,sa' no son necesarios.
 deallocate(ita,sa)
endif
deallocate(poi,podosi,potresi,ck)
end
	 
!-------------------------------------------------------------------------------------------------------------------------------------
!Subrutina CAJASASUBT
!A continuación se crea la caja As que aparece en la ecuación para flujo subterráneo (con esta ecuación la discretización es 
!siempre con elementos lineales, además todas las cajas irán en la misma posición). 
!'Nn' son las matrices elementales de esta caja. En 'vv' están los niveles freáticos solución de la iteración anterior de los que 
!se obtiene el espesor del que dependen estas cajas. Por ello se trata de una caja no lineal del sistema.
!También se da la condición de espesor mínimo en caso de que el espesor obtenido sea negativo.
!------------------------------------------------------------------------------------------------------------------------------------- 
subroutine cajasasubt(i,j,x,y,zp,vv,kix,kiy,ag)
use elemental
integer*4, dimension(:),allocatable::po
integer*4 i,j
real*8, dimension(:),allocatable::kxx,kyy,kxy
real*8 x(i),y(i),zp(i),Nn(3,3),Mpx(3),Mpy(3),Mh,P,kex,key,kexy,vv(i),kix(i),kiy(i),ag(i),a,b

allocate(po(i),kxx(i),kyy(i),kxy(i))
   
39     format(4/,A80)
40     format(6X,6(X,I5))

!Inizialización previa de variables:
!-----------------------------------
luno=(/0.5_8, 0.0_8, 0.5_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8/)		   
ldos=(/0.5_8, 0.5_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8/)
ja=(/1.0_8/6.0_8, 1.0_8/6.0_8, 1.0_8/6.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8/)
do u=1,i
po(u)=pos(u)
enddo
!Cálculo de las conductividades kxx, kyy y kxy a partir de las conductividades en las direcciones principales y el 
!ángulo de anisotropía.
do u=1,i
kxx(u)=kix(u)*(cos(ag(u))**2)+kiy(u)*(sin(ag(u))**2)
kyy(u)=kix(u)*(sin(ag(u))**2)+kiy(u)*(cos(ag(u))**2)
kxy(u)=(kix(u)-kiy(u))*sin(ag(u))*cos(ag(u))
enddo

!Cálculo de los valores para cada elemento
!-----------------------------------------
open(unit=3,file='C:\mallasub.txt',status='old')
read(3,39)ac
 do u=1,j
     read(3,40) n(1),n(4),n(2),n(5),n(3),n(6)
  
     xa=x(n(1))-x(n(3)) 
	 xb=x(n(2))-x(n(3)) 
	 ya=y(n(1))-y(n(3))
	 yb=y(n(2))-y(n(3))
	 jac=abs(xa*yb-xb*ya)

	 do uu=1,3
     do uuu=1,3
     Nn(uu,uuu)=0.0_8
     enddo
     enddo

	 do uu=1,3
	 Mpx(1)=(y(n(2))-y(n(3)))/jac
	 Mpx(2)=(y(n(3))-y(n(1)))/jac
	 Mpx(3)=(y(n(1))-y(n(2)))/jac
	 
	 Mpy(1)=(x(n(3))-x(n(2)))/jac
	 Mpy(2)=(x(n(1))-x(n(3)))/jac
	 Mpy(3)=(x(n(2))-x(n(1)))/jac  

	 !Se calculan valores en cada punto de integración del nivel freático, el sustrato impermeable y conductividades.
	 Mh=luno(uu)*vv(n(1))+ldos(uu)*vv(n(2))+(1.0_8-luno(uu)-ldos(uu))*vv(n(3))
	 P=luno(uu)*zp(n(1))+ldos(uu)*zp(n(2))+(1.0_8-luno(uu)-ldos(uu))*zp(n(3))
	 kex=luno(uu)*kxx(n(1))+ldos(uu)*kxx(n(2))+(1.0_8-luno(uu)-ldos(uu))*kxx(n(3))
	 key=luno(uu)*kyy(n(1))+ldos(uu)*kyy(n(2))+(1.0_8-luno(uu)-ldos(uu))*kyy(n(3))
	 kexy=luno(uu)*kxy(n(1))+ldos(uu)*kxy(n(2))+(1.0_8-luno(uu)-ldos(uu))*kxy(n(3))
	 
	 !Condición de espesor mínimo. En la ecuación de agua subterránea al igual que las ecuaciones de aguas someras no se pueden introducir valores
	 !negativos de espesor. Aquí no se modifica el valor de nivel freático de la iteración anterior, se reemplaza el espesor en el punto de 
	 !integración (Mh-P) por 0.001 si se obtiene un valor negativo. Es equivalente a que se modifique la conductividad en estos elementos de forma
	 !que la transmisividad sea positiva.
	 if ((Mh-P).le.0.0) then
	 Mh=P+0.001_8
	 endif 
	 		
	  do ui=1,3
  	  do uj=1,3
	  !Posicion 1,1. Cajas Asx+Asy.
	  a=kex*(Mh-P)*Mpx(ui)+kexy*(Mh-P)*Mpy(ui)
	  b=kexy*(Mh-P)*Mpx(ui)+key*(Mh-P)*Mpy(ui)
	  Nn(uj,ui)=Nn(uj,ui)+(Mpx(uj)*a+Mpy(uj)*b)*ja(uu)*jac
	  enddo
	  enddo	   
    enddo
	call suma(Nn,n,0,0,i,3,3,po)
enddo
close(3)
deallocate(po,kxx,kyy,kxy)
end
 
!------------------------------------------------------------------------------------------------------------------------------
!Subrutina TIMESUBT
!A continuación se calculan coeficientes de las matrices elementales que forman la caja Ns, típicamente llamada matriz de masa.
!'Nn' son las matrices elementales de esta caja. Estas cajas aparecen al discretizar la ecuación para flujo subterráneo.
!Sólo entra aquí en caso de resolver el problema no estacionario.
!------------------------------------------------------------------------------------------------------------------------------
subroutine timesubt(i,j,x,y,nd,At)
use elemental
integer*4, dimension(:),allocatable::po
integer*4 i,j				  
real*8 x(i),y(i),Nn(3,3),Mp(3),nd(i),Np,At

allocate(po(i))

39     format(4/,A80)
40     format(6X,6(X,I5))
					  
!Inizialización previa de variables:
!-----------------------------------
luno=(/1.0_8/3.0_8, 0.6_8, 0.2_8, 0.2_8, 0.0_8, 0.0_8, 0.0_8/)		   
ldos=(/1.0_8/3.0_8, 0.2_8, 0.6_8, 0.2_8, 0.0_8, 0.0_8, 0.0_8/)
ja=(/-27.0_8/96.0_8, 25.0_8/96.0_8, 25.0_8/96.0_8, 25.0_8/96.0_8, 0.0_8, 0.0_8, 0.0_8/)
do u=1,i
po(u)=pos(u)
enddo

!Cálculo de los valores para cada elemento
!-----------------------------------------
open(unit=3,file='C:\mallasub.txt',status='old')
read(3,39)ac	 
 do u=1,j
     read(3,40) n(1),n(4),n(2),n(5),n(3),n(6)
  
	 xa=x(n(1))-x(n(3)) 
	 xb=x(n(2))-x(n(3)) 
	 ya=y(n(1))-y(n(3))
	 yb=y(n(2))-y(n(3))  	 
	 jac=abs(xa*yb-xb*ya)
	 	 
	 do uu=1,3
     do uuu=1,3
     Nn(uu,uuu)=0.0_8
     enddo
     enddo

	 do uu=1,4			   
	 Mp(1)=luno(uu)
	 Mp(2)=ldos(uu)
	 Mp(3)=1.0_8-luno(uu)-ldos(uu)
	 
	 !Se calculan valores en cada punto de integración de la porosidad efectiva.
	 Np=Mp(1)*nd(n(1))+Mp(2)*nd(n(2))+Mp(3)*nd(n(3))

	  do ui=1,3
  	  do uj=1,3
      !Posición 1,1. Caja Ns.
	  Nn(uj,ui)=Nn(uj,ui)+(Mp(uj)*Np*Mp(ui))*ja(uu)*jac
	  enddo
	  enddo	  	 
	 enddo
	 call suma(Nn/At,n,0,0,i,3,3,po)
 enddo
 close(3)
 deallocate(po)
end

!--------------------------------------------------------------------------------------------------------------------------
!Subrutina TIMESUBTCONC
!Es una alternativa a la subrutina timesubt. Se diferencia de la anterior en que se concentra la masa en la diagonal.	
!En otras palabras, se suman todos los términos de cada fila en el término de la diagonal.
!--------------------------------------------------------------------------------------------------------------------------
subroutine timesubtconc(i,j,x,y,nd,At)
use elemental
integer*4, dimension(:),allocatable::po
integer*4 i,j				  
real*8 x(i),y(i),Nn(3,3),Mp(3),nd(i),Np,At

allocate(po(i))

39     format(4/,A80)
40     format(6X,6(X,I5))
				
!Inizialización previa de variables:
!-----------------------------------
luno=(/0.5_8, 0.0_8, 0.5_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8/)		   
ldos=(/0.5_8, 0.5_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8/)
ja=(/1.0_8/6.0_8, 1.0_8/6.0_8, 1.0_8/6.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8/)
do u=1,i
po(u)=pos(u)
enddo

!Cálculo de los valores para cada elemento
!-----------------------------------------
open(unit=3,file='C:\mallasub.txt',status='old')
read(3,39)ac	 
 do u=1,j
     read(3,40) n(1),n(4),n(2),n(5),n(3),n(6)
  
	 xa=x(n(1))-x(n(3)) 
	 xb=x(n(2))-x(n(3)) 
	 ya=y(n(1))-y(n(3))
	 yb=y(n(2))-y(n(3))  	 
	 jac=abs(xa*yb-xb*ya)
	 	 
	 do uu=1,3
     do uuu=1,3
     Nn(uu,uuu)=0.0_8
     enddo
     enddo

	 do uu=1,3
	 Mp(1)=luno(uu)
	 Mp(2)=ldos(uu)
	 Mp(3)=1.0_8-luno(uu)-ldos(uu)
	 
	 Np=Mp(1)*nd(n(1))+Mp(2)*nd(n(2))+Mp(3)*nd(n(3))
 
  	  do uj=1,3
      !Posición 1,1. Caja Ns concentrada.
	  !La suma de las funciones de forma en cualquier punto del elemento será igual a la unidad. 
	  Nn(uj,uj)=Nn(uj,uj)+(Mp(uj)*Np)*ja(uu)*jac
	  enddo	  	 
	 enddo
	 call suma(Nn/At,n,0,0,i,3,3,po)
 enddo
 close(3)
 deallocate(po)
end

!---------------------------------------------------------------------------------------------------------------------------------------------
!Subrutinas para el cálculo de coeficientes de las integrales en el dominio (formulación Petrov-Galerkin).
!Método de estabilización SUPG/PSPG con grad-div para las ecuaciones para flujo superficial. Parámetros de estabilización fueron desarrollados 
!por Gelhard para elementos P2-P1 para mallas uniformes considerando las ecuaciones de Navier-Stokes 2D en función de la presión cinemática y 
!sin considerar las integrales de contorno de las tensiones viscosas. La estabilización no está programada para el método de Newton.  
!Para elementos cuadráticos 'Mixx' serán las derivadas segundas en x de las funciones de forma y 'Miyy' las derivadas segundas en y de las 
!funciones de forma.
!---------------------------------------------------------------------------------------------------------------------------------------------
!Subrutina TIMEASUPGNS
!Para las ecuaciones de Navier-Stokes 2D. No se aplica como peso el término convectivo (despreciable) ni el término temporal.
!Se utiliza como longitud elemental el diámetro de la circunferencia que inscribe al elemento y unos parámetros similares a los de 
!Gelhard (aquí las ecuaciones están escritas en función del calado).
!La estabilización contiene en las dos primeras ecuaciones una ec. dinámica + ec. dinámica (con término supg)+ ec. cont (con término lsic)
!y como última ecuación la ec. continuidad + ec. dinámica (con peso pspg).
!Aquí se calculan las matrices elementales de los términos de masa estabilizados. Se utilizarían al resolver de forma transitoria (en el 
!esquema semi-implícito siempre) las ecuaciones, pero de momento no se utilizan.
!---------------------------------------------------------------------------------------------------------------------------------------------
subroutine timeasupgns (i,j,x,y,vv,At)
use elemental
integer*4, dimension(:),allocatable::poi,podosi
integer*4 i,j
real*8 x(i),y(i),Mi(6),Mix(6),Miy(6),Mpx(3),Mpy(3),Nt(6,6),Ntpx(3,6),Ntpy(3,6)	 
real*8 At,lsic,Mu,Mv,lunoq(13),ldosq(13),jaq(13),vv(3*i),cx,ccx,h,Muu,Mvv,a,b,c 

allocate(poi(3*i),podosi(3*i))

39     format(4/,A80)
40     format(6X,6(X,I5))
  
!Inizialización previa de variables:
!-----------------------------------
lunoq=(/0.0651301029002_8, 0.8697397941956_8, 0.0651301029002_8, 0.3128654960049_8, 0.6384441885698_8, 0.0486903154253_8, 	 &
0.6384441885698_8, 0.3128654960049_8, 0.0486903154253_8, 0.2603459660790_8, 0.4793080678419_8, 0.2603459660790_8, 1.0_8/3.0_8/)		   
ldosq=(/0.0651301029002_8, 0.0651301029002_8, 0.8697397941956_8, 0.0486903154253_8, 0.3128654960049_8, 0.6384441885698_8, 	 &
0.0486903154253_8, 0.6384441885698_8, 0.3128654960049_8, 0.2603459660790_8, 0.2603459660790_8, 0.4793080678419_8, 1.0_8/3.0_8/)
jaq=(/0.0266736178004_8, 0.0266736178004_8, 0.0266736178004_8, 0.03855688044515_8, 0.03855688044515_8, 0.03855688044515_8, 	 &
0.03855688044515_8, 0.03855688044515_8, 0.03855688044515_8, 0.0878076287166_8, 0.0878076287166_8, 0.0878076287166_8, 		 &
-0.0747850222335_8/)
do u=1,3*i
poi(u)=posi(u)
podosi(u)=posdosi(u)
enddo
lsic=1.0_8
!Introduciendo a=9.81 se consideraría la gravedad en el peso PSPG tal y como debería de ser si se utiliza el operador de la ecuación
!diferencial como peso. Introduciendo a=0.0 se consideraría la estabilización SUPG con grad-div. En cualquiera de los dos casos se 
!obtendrían peores soluciones.
!Se puede entender que se aplica la gravedad en el peso y un tercer parámetro ca=(h**2.0_8)/9.81_8 en vez de ccx donde se use 'a'
!de forma que los parámetros para SUPG y PSPG son diferentes.
a=1.0_8

!Cálculo de los valores para cada elemento
!-----------------------------------------
open(unit=1,file='C:\malla.txt',status='old')
read(1,39)ac
  do u=1,j
     read(1,40) n(1),n(4),n(2),n(5),n(3),n(6)

	 xa=x(n(1))-x(n(3)) 
	 xb=x(n(2))-x(n(3)) 
	 ya=y(n(1))-y(n(3))
	 yb=y(n(2))-y(n(3))  
	 jac=abs(xa*yb-xb*ya)
	 	 
	 do uu=1,6
     do uuu=1,6
	 Nt(uu,uuu)=0.0_8
     enddo
     enddo
	 do uu=1,3
     do uuu=1,6
	 Ntpx(uu,uuu)=0.0_8
	 Ntpy(uu,uuu)=0.0_8
     enddo
     enddo
     
	 Mpx(1)=(y(n(2))-y(n(3)))/jac
	 Mpx(2)=(y(n(3))-y(n(1)))/jac
	 Mpx(3)=(y(n(1))-y(n(2)))/jac
	 
	 Mpy(1)=(x(n(3))-x(n(2)))/jac
	 Mpy(2)=(x(n(1))-x(n(3)))/jac
	 Mpy(3)=(x(n(2))-x(n(1)))/jac
	 
	 !Valores medios de velocidad en el elemento (valores en el centro del elemento).
	 Muu=-(vv(n(1))+vv(n(2))+vv(n(3)))/9.0_8+(vv(n(4))+vv(n(5))+vv(n(6)))*4.0_8/9.0_8   
	 Mvv=-(vv(n(1)+i)+vv(n(2)+i)+vv(n(3)+i))/9.0_8+(vv(n(4)+i)+vv(n(5)+i)+vv(n(6)+i))*4.0_8/9.0_8

	 !Cálculo del diametro de la circunferencia inscrita:
	 a=sqrt((x(n(1))-x(n(2)))**2.0_8+(y(n(1))-y(n(2)))**2.0_8)
	 b=sqrt((x(n(2))-x(n(3)))**2.0_8+(y(n(2))-y(n(3)))**2.0_8)
	 c=sqrt((x(n(3))-x(n(1)))**2.0_8+(y(n(3))-y(n(1)))**2.0_8)
	 h=2.0_8*jac/(a+b+c)

	 !Parámetros de estabilización (del tipo de los de Tobias Geldar) para elementos que cumplen la condición LBB. Se podrían escribir 
	 !como cx=9.81_8, ccx=(h**2.0_8). Mismos parámetros para SUPG y PSPG (en ccx).
	 cx=9.81_8*lsic
	 ccx=9.81_8*(h**2.0_8)/cx

	 !Integración gauss de quinto orden con trece puntos de integración. Polinomios de grado superior a grado 7 quedarían 
	 !integrados con error.
	 do uu=1,13
	 
	 Mi(1)=lunoq(uu)*(2.0_8*lunoq(uu)-1.0_8)
	 Mi(2)=ldosq(uu)*(2.0_8*ldosq(uu)-1.0_8)
	 Mi(3)=(1.0_8-lunoq(uu)-ldosq(uu))*(1.0_8-2.0_8*lunoq(uu)-2.0_8*ldosq(uu))
	 Mi(4)=4.0_8*lunoq(uu)*ldosq(uu)
	 Mi(5)=4.0_8*ldosq(uu)*(1.0_8-lunoq(uu)-ldosq(uu))
	 Mi(6)=4.0_8*lunoq(uu)*(1.0_8-lunoq(uu)-ldosq(uu))

	 Mix(1)=((4.0_8*lunoq(uu)-1.0_8)*yb)/jac
	 Mix(2)=((1.0_8-4.0_8*ldosq(uu))*ya)/jac
	 Mix(3)=((4.0_8*lunoq(uu)+4.0_8*ldosq(uu)-3.0_8)*(yb-ya))/jac
	 Mix(4)=(4.0_8*ldosq(uu)*yb-4.0_8*lunoq(uu)*ya)/jac
	 Mix(5)=((8.0_8*ldosq(uu)+4.0_8*lunoq(uu)-4.0_8)*ya-4.0_8*ldosq(uu)*yb)/jac
	 Mix(6)=((4.0_8-4.0_8*ldosq(uu)-8.0_8*lunoq(uu))*yb+4.0_8*lunoq(uu)*ya)/jac
	
	 Miy(1)=((1.0_8-4.0_8*lunoq(uu))*xb)/jac
	 Miy(2)=((4.0_8*ldosq(uu)-1.0_8)*xa)/jac
	 Miy(3)=((4.0_8*lunoq(uu)+4.0_8*ldosq(uu)-3.0_8)*(xa-xb))/jac
	 Miy(4)=(4.0_8*lunoq(uu)*xa-4.0_8*ldosq(uu)*xb)/jac
	 Miy(5)=(4.0_8*ldosq(uu)*xb+(4.0_8-4.0_8*lunoq(uu)-8.0_8*ldosq(uu))*xa)/jac
	 Miy(6)=((8.0_8*lunoq(uu)+4.0_8*ldosq(uu)-4.0_8)*xb-4.0_8*lunoq(uu)*xa)/jac
 
	 Mu=dot_product(vv(n),Mi)   
	 Mv=dot_product(vv(n+i),Mi)

	  do ui=1,6
  	  do uj=1,6 
	  !Posición 1,1 y 2,2. Cajas de las matrices de masa M afectadas por el peso SUPG. 
	  Nt(uj,ui)=Nt(uj,ui)+ccx*(Mu*Mix(uj)+Mv*Miy(uj))*Mi(ui)*jaq(uu)*jac
	  enddo								  
	  enddo							 
	  
	  do ui=1,6	
	  do uj=1,3
	  !Posicion 3,3. Cajas de las matrices de masa M afectadas por el peso PSPG.	 	  	 	  
	  Ntpx(uj,ui)=Ntpx(uj,ui)+ccx*Mpx(uj)*Mi(ui)*jaq(uu)*jac
	  Ntpy(uj,ui)=Ntpy(uj,ui)+ccx*Mpy(uj)*Mi(ui)*jaq(uu)*jac 
	  enddo
	  enddo	    
	 enddo	 
	 call suma(Nt/At,n,0,0,3*i,6,6,poi)
	 call suma(Nt/At,n,i,i,3*i,6,6,podosi)	 
	 !Sin las siguientes matrices sólo se tendría la estabilización SUPG con grad-div.
	 call suma(a*Ntpx/At,n,2*i,0,3*i,3,6,poi)
	 call suma(a*Ntpy/At,n,2*i,i,3*i,3,6,podosi)
 enddo
 close(1)
 deallocate(poi,podosi)
end

!-----------------------------------------------------------------------------------------------------------------------
!Subrutina CAJASUPGNS
!Para las ecuaciones de Navier-Stokes 2D. Aquí se calculan las matrices elementales del resto de términos estabilizados.
!-----------------------------------------------------------------------------------------------------------------------
subroutine cajasupgns(i,j,x,y,vv,nu,del)
use elemental	  
integer*4, dimension(:),allocatable::poi,podosi,potresi
integer*4 i,j
real*8 x(i),y(i),Nn(6,6),Nm(6,6),Nk(6,3),Np(6,3),Mi(6),Mix(6),Miy(6),Mixx(6),Miyy(6),Mp(3),Mpx(3),Mpy(3)
real*8 Nh(3,6),Nf(3,6),Nl(3,3),Nj(3,6),Nw(3,6),Nxx(6,6),Nxy(6,6),Nyx(6,6),Nyy(6,6),lsic
real*8 Mu,Mv,lunoq(13),ldosq(13),jaq(13),vv(3*i),cx,ccx,nu,h,Muu,Mvv,a,b,c,del 

allocate(poi(3*i),podosi(3*i),potresi(3*i))

39     format(4/,A80)
40     format(6X,6(X,I5))

!Inizialización previa de variables:
!-----------------------------------
lunoq=(/0.0651301029002_8, 0.8697397941956_8, 0.0651301029002_8, 0.3128654960049_8, 0.6384441885698_8, 0.0486903154253_8, 	 &
0.6384441885698_8, 0.3128654960049_8, 0.0486903154253_8, 0.2603459660790_8, 0.4793080678419_8, 0.2603459660790_8, 1.0_8/3.0_8/)		   
ldosq=(/0.0651301029002_8, 0.0651301029002_8, 0.8697397941956_8, 0.0486903154253_8, 0.3128654960049_8, 0.6384441885698_8, 	 &
0.0486903154253_8, 0.6384441885698_8, 0.3128654960049_8, 0.2603459660790_8, 0.2603459660790_8, 0.4793080678419_8, 1.0_8/3.0_8/)
jaq=(/0.0266736178004_8, 0.0266736178004_8, 0.0266736178004_8, 0.03855688044515_8, 0.03855688044515_8, 0.03855688044515_8, 	 &
0.03855688044515_8, 0.03855688044515_8, 0.03855688044515_8, 0.0878076287166_8, 0.0878076287166_8, 0.0878076287166_8, 		 &
-0.0747850222335_8/)
do u=1,3*i
poi(u)=posi(u)
podosi(u)=posdosi(u)
potresi(u)=postresi(u)
enddo
lsic=1.0_8
a=1.0_8
  
!Cálculo de los valores para cada elemento
!-----------------------------------------
open(unit=1,file='C:\malla.txt',status='old')
read(1,39)ac
  do u=1,j
     read(1,40) n(1),n(4),n(2),n(5),n(3),n(6)

	 xa=x(n(1))-x(n(3)) 
	 xb=x(n(2))-x(n(3)) 
	 ya=y(n(1))-y(n(3))
	 yb=y(n(2))-y(n(3))  
	 jac=abs(xa*yb-xb*ya)
	 	 
	 do uu=1,6
     do uuu=1,6
     Nn(uu,uuu)=0.0_8
	 Nm(uu,uuu)=0.0_8	 
	 Nxx(uu,uuu)=0.0_8
	 Nxy(uu,uuu)=0.0_8
	 Nyx(uu,uuu)=0.0_8
	 Nyy(uu,uuu)=0.0_8
     enddo
     do uuu=1,3
	 Nk(uu,uuu)=0.0_8
	 Np(uu,uuu)=0.0_8	 
     enddo
     enddo
	 do uu=1,3
     do uuu=1,6
	 Nh(uu,uuu)=0.0_8
	 Nf(uu,uuu)=0.0_8	 
     enddo
     enddo
     
	 !Se calculan antes estos valores ya que son contantes para cada elemento.
	 Mixx(1)=4.0_8*(yb**2.0_8)/(jac**2.0_8)
	 Mixx(2)=4.0_8*(ya**2.0_8)/(jac**2.0_8)
	 Mixx(3)=4.0_8*((yb-ya)**2.0_8)/(jac**2.0_8)
	 Mixx(4)=-8.0_8*ya*yb/(jac**2.0_8)
	 Mixx(5)=8.0_8*(ya*yb-ya**2.0_8)/(jac**2.0_8)
	 Mixx(6)=8.0_8*(ya*yb-yb**2.0_8)/(jac**2.0_8)

	 Miyy(1)=4.0_8*(xb**2.0_8)/(jac**2.0_8)
	 Miyy(2)=4.0_8*(xa**2.0_8)/(jac**2.0_8)
	 Miyy(3)=4.0_8*((xa-xb)**2.0_8)/(jac**2.0_8)
	 Miyy(4)=-8.0_8*xa*xb/(jac**2.0_8)
	 Miyy(5)=8.0_8*(xa*xb-xa**2.0_8)/(jac**2.0_8)
	 Miyy(6)=8.0_8*(xa*xb-xb**2.0_8)/(jac**2.0_8)

	 Mpx(1)=(y(n(2))-y(n(3)))/jac
	 Mpx(2)=(y(n(3))-y(n(1)))/jac
	 Mpx(3)=(y(n(1))-y(n(2)))/jac
	 
	 Mpy(1)=(x(n(3))-x(n(2)))/jac
	 Mpy(2)=(x(n(1))-x(n(3)))/jac
	 Mpy(3)=(x(n(2))-x(n(1)))/jac
	 
	 Muu=-(vv(n(1))+vv(n(2))+vv(n(3)))/9.0_8+(vv(n(4))+vv(n(5))+vv(n(6)))*4.0_8/9.0_8   
	 Mvv=-(vv(n(1)+i)+vv(n(2)+i)+vv(n(3)+i))/9.0_8+(vv(n(4)+i)+vv(n(5)+i)+vv(n(6)+i))*4.0_8/9.0_8
	 
	 a=sqrt((x(n(1))-x(n(2)))**2.0_8+(y(n(1))-y(n(2)))**2.0_8)
	 b=sqrt((x(n(2))-x(n(3)))**2.0_8+(y(n(2))-y(n(3)))**2.0_8)
	 c=sqrt((x(n(3))-x(n(1)))**2.0_8+(y(n(3))-y(n(1)))**2.0_8)
	 h=2.0_8*jac/(a+b+c)

	 cx=9.81_8*lsic
	 ccx=9.81_8*(h**2.0_8)/cx
	 
	 do uu=1,13
     
	 Mi(1)=lunoq(uu)*(2.0_8*lunoq(uu)-1.0_8)
	 Mi(2)=ldosq(uu)*(2.0_8*ldosq(uu)-1.0_8)
	 Mi(3)=(1.0_8-lunoq(uu)-ldosq(uu))*(1.0_8-2.0_8*lunoq(uu)-2.0_8*ldosq(uu))
	 Mi(4)=4.0_8*lunoq(uu)*ldosq(uu)
	 Mi(5)=4.0_8*ldosq(uu)*(1.0_8-lunoq(uu)-ldosq(uu))
	 Mi(6)=4.0_8*lunoq(uu)*(1.0_8-lunoq(uu)-ldosq(uu))

	 Mix(1)=((4.0_8*lunoq(uu)-1.0_8)*yb)/jac
	 Mix(2)=((1.0_8-4.0_8*ldosq(uu))*ya)/jac
	 Mix(3)=((4.0_8*lunoq(uu)+4.0_8*ldosq(uu)-3.0_8)*(yb-ya))/jac
	 Mix(4)=(4.0_8*ldosq(uu)*yb-4.0_8*lunoq(uu)*ya)/jac
	 Mix(5)=((8.0_8*ldosq(uu)+4.0_8*lunoq(uu)-4.0_8)*ya-4.0_8*ldosq(uu)*yb)/jac
	 Mix(6)=((4.0_8-4.0_8*ldosq(uu)-8.0_8*lunoq(uu))*yb+4.0_8*lunoq(uu)*ya)/jac
	
	 Miy(1)=((1.0_8-4.0_8*lunoq(uu))*xb)/jac
	 Miy(2)=((4.0_8*ldosq(uu)-1.0_8)*xa)/jac
	 Miy(3)=((4.0_8*lunoq(uu)+4.0_8*ldosq(uu)-3.0_8)*(xa-xb))/jac
	 Miy(4)=(4.0_8*lunoq(uu)*xa-4.0_8*ldosq(uu)*xb)/jac
	 Miy(5)=(4.0_8*ldosq(uu)*xb+(4.0_8-4.0_8*lunoq(uu)-8.0_8*ldosq(uu))*xa)/jac
	 Miy(6)=((8.0_8*lunoq(uu)+4.0_8*ldosq(uu)-4.0_8)*xb-4.0_8*lunoq(uu)*xa)/jac
 
	 Mp(1)=lunoq(uu)
	 Mp(2)=ldosq(uu)
	 Mp(3)=1.0_8-lunoq(uu)-ldosq(uu)

	 Mu=dot_product(vv(n),Mi)   
	 Mv=dot_product(vv(n+i),Mi)

	  do ui=1,6
  	  do uj=1,6 
	  !Posición 1,1 y 2,2. 
	  !Cajas no lineales C con peso SUPG. 
	  Nn(uj,ui)=Nn(uj,ui)+ccx*(Mu*Mix(uj)+Mv*Miy(uj))*(Mix(ui)*Mu+Miy(ui)*Mv)*jaq(uu)*jac	  
	  !Caja A con peso supg (sin forma débil)
	  Nm(uj,ui)=Nm(uj,ui)+ccx*(Mu*Mix(uj)+Mv*Miy(uj))*(Miyy(ui)+Mixx(ui))*jaq(uu)*jac	  
	  !Cajas correspondientes a la ecuación de continuidad estabilizada con pesos grad-div (se usa el parámetro 'cx').
	  Nxx(uj,ui)=Nxx(uj,ui)+cx*Mix(uj)*Mix(ui)*jaq(uu)*jac
	  Nxy(uj,ui)=Nxy(uj,ui)+cx*Mix(uj)*Miy(ui)*jaq(uu)*jac	  
	  Nyx(uj,ui)=Nyx(uj,ui)+cx*Miy(uj)*Mix(ui)*jaq(uu)*jac
	  Nyy(uj,ui)=Nyy(uj,ui)+cx*Miy(uj)*Miy(ui)*jaq(uu)*jac	  
	  enddo								  
	  enddo	
	  
	  do uj=1,6
	  do ui=1,3
	  !Posicion 1,3. Caja Bx con peso SUPG.
	  !Sin forma débil	 	  
	  Nk(uj,ui)=Nk(uj,ui)+ccx*(Mu*Mix(uj)+Mv*Miy(uj))*Mpx(ui)*jaq(uu)*jac	  
	  !Posicion 2,3. Caja By con peso SUPG.
	  !Sin forma débil
	  Np(uj,ui)=Np(uj,ui)+ccx*(Mu*Mix(uj)+Mv*Miy(uj))*Mpy(ui)*jaq(uu)*jac	 
	  enddo
	  enddo
	  
	  do ui=1,6	
	  do uj=1,3
	  !Posicion 3,1. Caja C con peso SUPG.	 	  
	  Nh(uj,ui)=Nh(uj,ui)+ccx*Mpx(uj)*(Mix(ui)*Mu+Miy(ui)*Mv)*jaq(uu)*jac
	  !Posicion 3,2. Caja C con peso SUPG.	 	  
	  Nf(uj,ui)=Nf(uj,ui)+ccx*Mpy(uj)*(Mix(ui)*Mu+Miy(ui)*Mv)*jaq(uu)*jac	  	  
	  enddo
	  enddo	     
	 enddo
	 
	 !Las siguientes cajas (con pesos PSPG) no se integran ya que contienen coeficientes contantes para cada elemento.	 
	 do ui=1,6	
	  do uj=1,3
	  !Posicion 3,1. Caja A con peso PSPG.	 	  
	  Nj(uj,ui)=ccx*Mpx(uj)*(Miyy(ui)+Mixx(ui))*jac/2.0_8
	  !Posicion 3,2. Caja A con peso PSPG.	 	  
	  Nw(uj,ui)=ccx*Mpy(uj)*(Miyy(ui)+Mixx(ui))*jac/2.0_8	  	  
	  enddo
	 enddo
	 do ui=1,3	
	  do uj=1,3
	  !Posicion 3,3. Caja Bx+By con peso PSPG.	 	  
	  Nl(uj,ui)=ccx*(Mpx(uj)*Mpx(ui)+Mpy(uj)*Mpy(ui))*jac/2.0_8	 	   
	  enddo
	 enddo
	 call suma((Nn-nu*Nm+Nxx)/del,n,0,0,3*i,6,6,poi)
	 call suma(Nxy/del,n,0,i,3*i,6,6,podosi)
	 call suma(Nyx/del,n,i,0,3*i,6,6,poi)
	 call suma((Nn-nu*Nm+Nyy)/del,n,i,i,3*i,6,6,podosi)
	 call suma(9.81_8*Nk/del,n,0,2*i,3*i,6,3,potresi)
	 call suma(9.81_8*Np/del,n,i,2*i,3*i,6,3,potresi)	 
	 !Sin las siguientes matrices sólo se tendría la estabilización SUPG con grad-div.
	 call suma((a*Nh-a*nu*Nj)/del,n,2*i,0,3*i,3,6,poi)
	 call suma((a*Nf-a*nu*Nw)/del,n,2*i,i,3*i,3,6,podosi)
	 call suma(a*9.81_8*Nl/del,n,2*i,2*i,3*i,3,3,potresi)
 enddo
 close(1)
 deallocate(poi,podosi,potresi)
end

!------------------------------------------------------------------------------------------------------------------------------------------  
!Subrutina TIMEASUPGAS
!Para las ecuaciones de aguas someras. No se aplica como peso el término convectivo (despreciable) ni el término de maning ni el 
!término temporal. De momento la estabilización no está aplicada al término de lluvia (iría en término independiente de ec dinámicas 
!con peso cx).
!Se utiliza como longitud elemental el diámetro de la circunferencia que inscribe al elemento y unos parámetros similares a los de la 
!estabilización de las ecuaciones de Navier-Stokes 2D.
!La estabilización contiene en las dos primeras ecuaciones una ec. dinámica + ec. dinámica (con término supg)+ ec. cont (con término lsic)
!y como última ecuación la ec. continuidad + ec. dinámica (con peso pspg).
!Aquí se calculan las matrices elementales de los términos de masa estabilizados. Se utilizarían al resolver de forma transitoria (en el 
!esquema semi-implícito siempre) las ecuaciones, pero de momento no se utilizan.
!-----------------------------------------------------------------------------------------------------------------------------------------
subroutine timeasupgas (i,j,x,y,vv,At)
use elemental
integer*4, dimension(:),allocatable::poi,podosi,potresi
integer*4 i,j
real*8 x(i),y(i),Mi(6),Mix(6),Miy(6),Mp(3),Mpx(3),Mpy(3),Nt(6,6),Ntpx(3,6),Ntpy(3,6)
real*8 At,lsic,Mu,Mv,lunoq(13),ldosq(13),jaq(13),vv(3*i),cx,ccx,h,Muu,Mvv,Mhu,a,b,c,hxy,Mhx,Mhy,Ntx(6,3),Nty(6,3)

allocate(poi(3*i),podosi(3*i),potresi(3*i))

39     format(4/,A80)
40     format(6X,6(X,I5))

!Inizialización previa de variables:
!-----------------------------------
lunoq=(/0.0651301029002_8, 0.8697397941956_8, 0.0651301029002_8, 0.3128654960049_8, 0.6384441885698_8, 0.0486903154253_8, 	 &
0.6384441885698_8, 0.3128654960049_8, 0.0486903154253_8, 0.2603459660790_8, 0.4793080678419_8, 0.2603459660790_8, 1.0_8/3.0_8/)		   
ldosq=(/0.0651301029002_8, 0.0651301029002_8, 0.8697397941956_8, 0.0486903154253_8, 0.3128654960049_8, 0.6384441885698_8, 	 &
0.0486903154253_8, 0.6384441885698_8, 0.3128654960049_8, 0.2603459660790_8, 0.2603459660790_8, 0.4793080678419_8, 1.0_8/3.0_8/)
jaq=(/0.0266736178004_8, 0.0266736178004_8, 0.0266736178004_8, 0.03855688044515_8, 0.03855688044515_8, 0.03855688044515_8, 	 &
0.03855688044515_8, 0.03855688044515_8, 0.03855688044515_8, 0.0878076287166_8, 0.0878076287166_8, 0.0878076287166_8, 		 &
-0.0747850222335_8/)
do u=1,3*i
poi(u)=posi(u)
podosi(u)=posdosi(u)
potresi(u)=postresi(u)
enddo
!Las consideraciones sobre la variable 'a' son idénticas a las explicadas en la estabilización de la ecuación de Navier-Stokes 2D.
lsic=1.0_8
a=1.0_8
  
!Cálculo de los valores para cada elemento
!-----------------------------------------
open(unit=1,file='C:\malla.txt',status='old')
read(1,39)ac
  do u=1,j
     read(1,40) n(1),n(4),n(2),n(5),n(3),n(6)

	 xa=x(n(1))-x(n(3)) 
	 xb=x(n(2))-x(n(3)) 
	 ya=y(n(1))-y(n(3))
	 yb=y(n(2))-y(n(3))  
	 jac=abs(xa*yb-xb*ya)
	 	 
	 do uu=1,6
     do uuu=1,6
	 Nt(uu,uuu)=0.0_8
     enddo
	 do uuu=1,3
	 Ntx(uu,uuu)=0.0_8
	 Nty(uu,uuu)=0.0_8
     enddo
     enddo
	 do uu=1,3
     do uuu=1,6
	 Ntpx(uu,uuu)=0.0_8
	 Ntpy(uu,uuu)=0.0_8
     enddo
     enddo
     
	 Mpx(1)=(y(n(2))-y(n(3)))/jac
	 Mpx(2)=(y(n(3))-y(n(1)))/jac
	 Mpx(3)=(y(n(1))-y(n(2)))/jac
	 
	 Mpy(1)=(x(n(3))-x(n(2)))/jac
	 Mpy(2)=(x(n(1))-x(n(3)))/jac
	 Mpy(3)=(x(n(2))-x(n(1)))/jac

	 !'Mhx, Mhy' son las derivadas del calado en dirección x y en dirección y evaluadas en cualquier punto de integración del elemento, 
	 !a partir de los valores de calado de la iteración anterior. Su valor es constante en todos los puntos.
	 Mhx=vv(2*i+n(1))*Mpx(1)+vv(2*i+n(2))*Mpx(2)+vv(2*i+n(3))*Mpx(3)
	 Mhy=vv(2*i+n(1))*Mpy(1)+vv(2*i+n(2))*Mpy(2)+vv(2*i+n(3))*Mpy(3)

	 !Valores medios de velocidad y calado en el elemento (valores en el centro del elemento).
	 Muu=-(vv(n(1))+vv(n(2))+vv(n(3)))/9.0_8+(vv(n(4))+vv(n(5))+vv(n(6)))*4.0_8/9.0_8   
	 Mvv=-(vv(n(1)+i)+vv(n(2)+i)+vv(n(3)+i))/9.0_8+(vv(n(4)+i)+vv(n(5)+i)+vv(n(6)+i))*4.0_8/9.0_8
	 Mhu=(vv(n(1)+2*i)+vv(n(2)+2*i)+vv(n(3)+2*i))/3.0_8

	 !Cálculo del diametro de la circunferencia inscrita:
	 a=sqrt((x(n(1))-x(n(2)))**2.0_8+(y(n(1))-y(n(2)))**2.0_8)
	 b=sqrt((x(n(2))-x(n(3)))**2.0_8+(y(n(2))-y(n(3)))**2.0_8)
	 c=sqrt((x(n(3))-x(n(1)))**2.0_8+(y(n(3))-y(n(1)))**2.0_8)
	 h=2.0_8*jac/(a+b+c)

	 !Parámetros de estabilización en función de 'Mhu' para elementos que cumplen la condición LBB.	Se podrían escribir 
	 !como cx=9.81_8/(Mhu**2.0_8), ccx=(h**2.0_8). Se utilizan los mismos parámetros para SUPG y PSPG (en ccx).  
	 cx=9.81_8*lsic/(Mhu**2.0_8)
	 ccx=9.81_8*(h**2.0_8)/(cx*(Mhu**2.0_8))
	 
	 do uu=1,13
     
	 Mi(1)=lunoq(uu)*(2.0_8*lunoq(uu)-1.0_8)
	 Mi(2)=ldosq(uu)*(2.0_8*ldosq(uu)-1.0_8)
	 Mi(3)=(1.0_8-lunoq(uu)-ldosq(uu))*(1.0_8-2.0_8*lunoq(uu)-2.0_8*ldosq(uu))
	 Mi(4)=4.0_8*lunoq(uu)*ldosq(uu)
	 Mi(5)=4.0_8*ldosq(uu)*(1.0_8-lunoq(uu)-ldosq(uu))
	 Mi(6)=4.0_8*lunoq(uu)*(1.0_8-lunoq(uu)-ldosq(uu))

	 Mix(1)=((4.0_8*lunoq(uu)-1.0_8)*yb)/jac
	 Mix(2)=((1.0_8-4.0_8*ldosq(uu))*ya)/jac
	 Mix(3)=((4.0_8*lunoq(uu)+4.0_8*ldosq(uu)-3.0_8)*(yb-ya))/jac
	 Mix(4)=(4.0_8*ldosq(uu)*yb-4.0_8*lunoq(uu)*ya)/jac
	 Mix(5)=((8.0_8*ldosq(uu)+4.0_8*lunoq(uu)-4.0_8)*ya-4.0_8*ldosq(uu)*yb)/jac
	 Mix(6)=((4.0_8-4.0_8*ldosq(uu)-8.0_8*lunoq(uu))*yb+4.0_8*lunoq(uu)*ya)/jac
	
	 Miy(1)=((1.0_8-4.0_8*lunoq(uu))*xb)/jac
	 Miy(2)=((4.0_8*ldosq(uu)-1.0_8)*xa)/jac
	 Miy(3)=((4.0_8*lunoq(uu)+4.0_8*ldosq(uu)-3.0_8)*(xa-xb))/jac
	 Miy(4)=(4.0_8*lunoq(uu)*xa-4.0_8*ldosq(uu)*xb)/jac
	 Miy(5)=(4.0_8*ldosq(uu)*xb+(4.0_8-4.0_8*lunoq(uu)-8.0_8*ldosq(uu))*xa)/jac
	 Miy(6)=((8.0_8*lunoq(uu)+4.0_8*ldosq(uu)-4.0_8)*xb-4.0_8*lunoq(uu)*xa)/jac
 
	 Mp(1)=lunoq(uu)
	 Mp(2)=ldosq(uu)
	 Mp(3)=1.0_8-lunoq(uu)-ldosq(uu)

	 Mu=dot_product(vv(n),Mi)   
	 Mv=dot_product(vv(n+i),Mi)
	 hxy=vv(2*i+n(1))*Mp(1)+vv(2*i+n(2))*Mp(2)+vv(2*i+n(3))*Mp(3)
	 
	  do ui=1,6
  	  do uj=1,6 
	  !Posición 1,1 y 2,2. Cajas de las matrices de masa M afectadas por el peso SUPG. 
	  Nt(uj,ui)=Nt(uj,ui)+ccx*(Mu*Mix(uj)+Mv*Miy(uj))*Mi(ui)*jaq(uu)*jac
	  enddo								  
	  enddo	
	  
	  do uj=1,6
	  do ui=1,3
	  !Cajas de las matrices de masa N afectadas por el peso grad-div. Estas no aparecen si se usan las ecuaciones de Navier-Stokes 2D.
	  !Posicion 1,3.
	  Ntx(uj,ui)=Ntx(uj,ui)+cx*(hxy*Mix(uj)+Mhx*Mi(uj))*Mp(ui)*jaq(uu)*jac	  
	  !Posicion 2,3.
	  Nty(uj,ui)=Nty(uj,ui)+cx*(hxy*Miy(uj)+Mhy*Mi(uj))*Mp(ui)*jaq(uu)*jac	 
	  enddo
	  enddo
  
	  do ui=1,6	
	  do uj=1,3
	  !Posicion 3,3. Cajas de las matrices de masa M afectadas por el peso PSPG.	 	  	 	  
	  Ntpx(uj,ui)=Ntpx(uj,ui)+ccx*Mpx(uj)*Mi(ui)*jaq(uu)*jac
	  Ntpy(uj,ui)=Ntpy(uj,ui)+ccx*Mpy(uj)*Mi(ui)*jaq(uu)*jac 
	  enddo
	  enddo	    
	 enddo	 
	 call suma(Nt/At,n,0,0,3*i,6,6,poi)
	 call suma(Nt/At,n,i,i,3*i,6,6,podosi)	 
	 !Sin las siguientes matrices sólo se tendría la estabilización SUPG con grad-div.
	 call suma(Ntx/At,n,0,2*i,3*i,6,3,potresi)  
	 call suma(Nty/At,n,i,2*i,3*i,6,3,potresi)  
	 call suma(a*Ntpx/At,n,2*i,0,3*i,3,6,poi)
	 call suma(a*Ntpy/At,n,2*i,i,3*i,3,6,podosi)
 enddo
 close(1)
deallocate(poi,podosi,potresi)
end

!-----------------------------------------------------------------------------------------------------------------------------------------
!Subrutina CAJASUPGAS
!Para las ecuaciones de aguas someras. Se calculan las matrices elementales de otros términos estabilizados.
!-----------------------------------------------------------------------------------------------------------------------------------------
subroutine cajasupgas(i,j,x,y,vv,nu,del)
use elemental
integer*4, dimension(:),allocatable::poi,podosi,potresi	  
integer*4 i,j
real*8 x(i),y(i),Nn(6,6),Nm(6,6),Nk(6,3),Np(6,3),Mi(6),Mix(6),Miy(6),Mixx(6),Miyy(6),Mp(3),Mpx(3),Mpy(3)
real*8 Nh(3,6),Nf(3,6),Nl(3,3),Nj(3,6),Nw(3,6),Nxx(6,6),Nxy(6,6),Nyx(6,6),Nyy(6,6),lsic
real*8 Mu,Mv,lunoq(13),ldosq(13),jaq(13),vv(3*i),cx,ccx,nu,h,Muu,Mvv,a,b,c,hxy,Mhx,Mhy,del,Mhu 

allocate(poi(3*i),podosi(3*i),potresi(3*i))

39     format(4/,A80)
40     format(6X,6(X,I5))

!Inizialización previa de variables:
!-----------------------------------
lunoq=(/0.0651301029002_8, 0.8697397941956_8, 0.0651301029002_8, 0.3128654960049_8, 0.6384441885698_8, 0.0486903154253_8, 	 &
0.6384441885698_8, 0.3128654960049_8, 0.0486903154253_8, 0.2603459660790_8, 0.4793080678419_8, 0.2603459660790_8, 1.0_8/3.0_8/)		   
ldosq=(/0.0651301029002_8, 0.0651301029002_8, 0.8697397941956_8, 0.0486903154253_8, 0.3128654960049_8, 0.6384441885698_8, 	 &
0.0486903154253_8, 0.6384441885698_8, 0.3128654960049_8, 0.2603459660790_8, 0.2603459660790_8, 0.4793080678419_8, 1.0_8/3.0_8/)
jaq=(/0.0266736178004_8, 0.0266736178004_8, 0.0266736178004_8, 0.03855688044515_8, 0.03855688044515_8, 0.03855688044515_8, 	 &
0.03855688044515_8, 0.03855688044515_8, 0.03855688044515_8, 0.0878076287166_8, 0.0878076287166_8, 0.0878076287166_8, 		 &
-0.0747850222335_8/)
do u=1,3*i
poi(u)=posi(u)
podosi(u)=posdosi(u)
potresi(u)=postresi(u)
enddo
lsic=1.0_8
a=1.0_8
  
!Cálculo de los valores para cada elemento
!-----------------------------------------
open(unit=1,file='C:\malla.txt',status='old')
read(1,39)ac
  do u=1,j
     read(1,40) n(1),n(4),n(2),n(5),n(3),n(6)

	 xa=x(n(1))-x(n(3)) 
	 xb=x(n(2))-x(n(3)) 
	 ya=y(n(1))-y(n(3))
	 yb=y(n(2))-y(n(3))  
	 jac=abs(xa*yb-xb*ya)
	 	 
	 do uu=1,6
     do uuu=1,6
     Nn(uu,uuu)=0.0_8
	 Nm(uu,uuu)=0.0_8
	 Nxx(uu,uuu)=0.0_8
	 Nxy(uu,uuu)=0.0_8
	 Nyx(uu,uuu)=0.0_8
	 Nyy(uu,uuu)=0.0_8
     enddo
	 do uuu=1,3
	 Nk(uu,uuu)=0.0_8
	 Np(uu,uuu)=0.0_8	 
     enddo
     enddo
	 do uu=1,3
     do uuu=1,6
	 Nh(uu,uuu)=0.0_8
	 Nf(uu,uuu)=0.0_8
     enddo
     enddo

	 Mixx(1)=4.0_8*(yb**2.0_8)/(jac**2.0_8)
	 Mixx(2)=4.0_8*(ya**2.0_8)/(jac**2.0_8)
	 Mixx(3)=4.0_8*((yb-ya)**2.0_8)/(jac**2.0_8)
	 Mixx(4)=-8.0_8*ya*yb/(jac**2.0_8)
	 Mixx(5)=8.0_8*(ya*yb-ya**2.0_8)/(jac**2.0_8)
	 Mixx(6)=8.0_8*(ya*yb-yb**2.0_8)/(jac**2.0_8)

	 Miyy(1)=4.0_8*(xb**2.0_8)/(jac**2.0_8)
	 Miyy(2)=4.0_8*(xa**2.0_8)/(jac**2.0_8)
	 Miyy(3)=4.0_8*((xa-xb)**2.0_8)/(jac**2.0_8)
	 Miyy(4)=-8.0_8*xa*xb/(jac**2.0_8)
	 Miyy(5)=8.0_8*(xa*xb-xa**2.0_8)/(jac**2.0_8)
	 Miyy(6)=8.0_8*(xa*xb-xb**2.0_8)/(jac**2.0_8)

	 Mpx(1)=(y(n(2))-y(n(3)))/jac
	 Mpx(2)=(y(n(3))-y(n(1)))/jac
	 Mpx(3)=(y(n(1))-y(n(2)))/jac
	 
	 Mpy(1)=(x(n(3))-x(n(2)))/jac
	 Mpy(2)=(x(n(1))-x(n(3)))/jac
	 Mpy(3)=(x(n(2))-x(n(1)))/jac

	 Mhx=vv(2*i+n(1))*Mpx(1)+vv(2*i+n(2))*Mpx(2)+vv(2*i+n(3))*Mpx(3)
	 Mhy=vv(2*i+n(1))*Mpy(1)+vv(2*i+n(2))*Mpy(2)+vv(2*i+n(3))*Mpy(3)

	 Muu=-(vv(n(1))+vv(n(2))+vv(n(3)))/9.0_8+(vv(n(4))+vv(n(5))+vv(n(6)))*4.0_8/9.0_8   
	 Mvv=-(vv(n(1)+i)+vv(n(2)+i)+vv(n(3)+i))/9.0_8+(vv(n(4)+i)+vv(n(5)+i)+vv(n(6)+i))*4.0_8/9.0_8
	 Mhu=(vv(n(1)+2*i)+vv(n(2)+2*i)+vv(n(3)+2*i))/3.0_8

	 a=sqrt((x(n(1))-x(n(2)))**2.0_8+(y(n(1))-y(n(2)))**2.0_8)
	 b=sqrt((x(n(2))-x(n(3)))**2.0_8+(y(n(2))-y(n(3)))**2.0_8)
	 c=sqrt((x(n(3))-x(n(1)))**2.0_8+(y(n(3))-y(n(1)))**2.0_8)
	 h=2.0_8*jac/(a+b+c)
 
	 cx=9.81_8*lsic/(Mhu**2.0_8)
	 ccx=9.81_8*(h**2.0_8)/(cx*(Mhu**2.0_8))  
	 
	 do uu=1,13
     
	 Mi(1)=lunoq(uu)*(2.0_8*lunoq(uu)-1.0_8)
	 Mi(2)=ldosq(uu)*(2.0_8*ldosq(uu)-1.0_8)
	 Mi(3)=(1.0_8-lunoq(uu)-ldosq(uu))*(1.0_8-2.0_8*lunoq(uu)-2.0_8*ldosq(uu))
	 Mi(4)=4.0_8*lunoq(uu)*ldosq(uu)
	 Mi(5)=4.0_8*ldosq(uu)*(1.0_8-lunoq(uu)-ldosq(uu))
	 Mi(6)=4.0_8*lunoq(uu)*(1.0_8-lunoq(uu)-ldosq(uu))

	 Mix(1)=((4.0_8*lunoq(uu)-1.0_8)*yb)/jac
	 Mix(2)=((1.0_8-4.0_8*ldosq(uu))*ya)/jac
	 Mix(3)=((4.0_8*lunoq(uu)+4.0_8*ldosq(uu)-3.0_8)*(yb-ya))/jac
	 Mix(4)=(4.0_8*ldosq(uu)*yb-4.0_8*lunoq(uu)*ya)/jac
	 Mix(5)=((8.0_8*ldosq(uu)+4.0_8*lunoq(uu)-4.0_8)*ya-4.0_8*ldosq(uu)*yb)/jac
	 Mix(6)=((4.0_8-4.0_8*ldosq(uu)-8.0_8*lunoq(uu))*yb+4.0_8*lunoq(uu)*ya)/jac
	
	 Miy(1)=((1.0_8-4.0_8*lunoq(uu))*xb)/jac
	 Miy(2)=((4.0_8*ldosq(uu)-1.0_8)*xa)/jac
	 Miy(3)=((4.0_8*lunoq(uu)+4.0_8*ldosq(uu)-3.0_8)*(xa-xb))/jac
	 Miy(4)=(4.0_8*lunoq(uu)*xa-4.0_8*ldosq(uu)*xb)/jac
	 Miy(5)=(4.0_8*ldosq(uu)*xb+(4.0_8-4.0_8*lunoq(uu)-8.0_8*ldosq(uu))*xa)/jac
	 Miy(6)=((8.0_8*lunoq(uu)+4.0_8*ldosq(uu)-4.0_8)*xb-4.0_8*lunoq(uu)*xa)/jac
 
	 Mp(1)=lunoq(uu)
	 Mp(2)=ldosq(uu)
	 Mp(3)=1.0_8-lunoq(uu)-ldosq(uu)

	 Mu=dot_product(vv(n),Mi)   
	 Mv=dot_product(vv(n+i),Mi)
	 hxy=vv(2*i+n(1))*Mp(1)+vv(2*i+n(2))*Mp(2)+vv(2*i+n(3))*Mp(3)
	 
	  do ui=1,6
  	  do uj=1,6 
	  !Posición 1,1 y 2,2. Cajas no lineales C con peso SUPG. 
	  Nn(uj,ui)=Nn(uj,ui)+ccx*(Mu*Mix(uj)+Mv*Miy(uj))*(Mix(ui)*Mu+Miy(ui)*Mv)*jaq(uu)*jac	  
	  !Caja A con peso SUPG (sin forma débil).
	  Nm(uj,ui)=Nm(uj,ui)+ccx*(Mu*Mix(uj)+Mv*Miy(uj))*(Miyy(ui)+Mixx(ui))*jaq(uu)*jac	  
	  !Cajas correspondientes a la ecuación de continuidad estabilizada con pesos grad-div (se usa el parámetro 'cx'). Obviamente, son 
	  !diferentes de las usadas para las ecuaciones de Navier-Stokes 2D.
	  Nxx(uj,ui)=Nxx(uj,ui)+cx*(hxy*Mix(uj)+Mhx*Mi(uj))*(hxy*Mix(ui)+Mhx*Mi(ui))*jaq(uu)*jac	 
	  Nxy(uj,ui)=Nxy(uj,ui)+cx*(hxy*Mix(uj)+Mhx*Mi(uj))*(hxy*Miy(ui)+Mhy*Mi(ui))*jaq(uu)*jac	  
	  Nyx(uj,ui)=Nyx(uj,ui)+cx*(hxy*Miy(uj)+Mhy*Mi(uj))*(hxy*Mix(ui)+Mhx*Mi(ui))*jaq(uu)*jac	 
	  Nyy(uj,ui)=Nyy(uj,ui)+cx*(hxy*Miy(uj)+Mhy*Mi(uj))*(hxy*Miy(ui)+Mhy*Mi(ui))*jaq(uu)*jac	 
	  enddo								  
	  enddo	
	  
	  do uj=1,6
	  do ui=1,3
	  !Posicion 1,3. Caja Bx con peso SUPG.
	  !Sin forma débil	 	  
	  Nk(uj,ui)=Nk(uj,ui)+ccx*(Mu*Mix(uj)+Mv*Miy(uj))*Mpx(ui)*jaq(uu)*jac	  
	  !Posicion 2,3. Caja By con peso SUPG.
	  !Sin forma débil
	  Np(uj,ui)=Np(uj,ui)+ccx*(Mu*Mix(uj)+Mv*Miy(uj))*Mpy(ui)*jaq(uu)*jac	 
	  enddo
	  enddo
	  
	  do ui=1,6	
	  do uj=1,3
	  !Posicion 3,1. Caja C con peso PSPG.	 	  
	  Nh(uj,ui)=Nh(uj,ui)+ccx*Mpx(uj)*(Mix(ui)*Mu+Miy(ui)*Mv)*jaq(uu)*jac
	  !Posicion 3,2. Caja C con peso PSPG.	 	  
	  Nf(uj,ui)=Nf(uj,ui)+ccx*Mpy(uj)*(Mix(ui)*Mu+Miy(ui)*Mv)*jaq(uu)*jac	  	  
	  enddo
	  enddo	   
	 enddo
	 
	 !Las siguientes cajas (con pesos PSPG) no se integran ya que contienen coeficientes contantes para cada elemento	 
	 do ui=1,6	
	  do uj=1,3
	  !Posicion 3,1. Caja A con peso PSPG.	 	  
	  Nj(uj,ui)=ccx*Mpx(uj)*(Miyy(ui)+Mixx(ui))*jac/2.0_8
	  !Posicion 3,2. Caja A con peso PSPG.	 	  
	  Nw(uj,ui)=ccx*Mpy(uj)*(Miyy(ui)+Mixx(ui))*jac/2.0_8	  	  
	  enddo
	 enddo
	 do ui=1,3	
	  do uj=1,3
	  !Posicion 3,3. Caja Bx+By con peso PSPG.	 	  
	  Nl(uj,ui)=ccx*(Mpx(uj)*Mpx(ui)+Mpy(uj)*Mpy(ui))*jac/2.0_8	 	   
	  enddo
	 enddo
	 call suma((Nn-nu*Nm+Nxx)/del,n,0,0,3*i,6,6,poi)
	 call suma(Nxy/del,n,0,i,3*i,6,6,podosi)
	 call suma(Nyx/del,n,i,0,3*i,6,6,poi)
	 call suma((Nn-nu*Nm+Nyy)/del,n,i,i,3*i,6,6,podosi)
	 call suma(9.81_8*Nk/del,n,0,2*i,3*i,6,3,potresi)
	 call suma(9.81_8*Np/del,n,i,2*i,3*i,6,3,potresi)	 
	 !Sin las siguientes matrices sólo se tendría la estabilización SUPG con grad-div.
	 call suma((a*Nh-a*nu*Nj)/del,n,2*i,0,3*i,3,6,poi)
	 call suma((a*Nf-a*nu*Nw)/del,n,2*i,i,3*i,3,6,podosi)
	 call suma(a*9.81_8*Nl/del,n,2*i,2*i,3*i,3,3,potresi)
 enddo
 close(1)
 deallocate(poi,podosi,potresi)
end

!-----------------------------------------------------------------------------------------------------------------------------------------
!Subrutina FSUPG
!Para las ecuaciones de aguas someras. En esta subrutina se calcula el término de fricción (de Manning) estabilizado.
!Se calcula para el caso yn='no' indicado en la subrutina f que evalúa las ecuaciones de aguas someras tal y como son.
!Se calculan vectores elementales que irán al término independiente del sistema.
!-----------------------------------------------------------------------------------------------------------------------------------------
subroutine fsupg(i,j,x,y,ma,vv,Nx,Ny,Nz) 
use elemental
integer*4 i,j	 
real*8 x(i),y(i),Nx(i),Ny(i),s,maning,ma(i),vv(3*i),Mi(6),Mu,Mv,Mui,Mvi,Mhi,dist,a,b,c
real*8 lsic,ccx,aa,Nz(i),Mix(6),Miy(6),Mpx(3),Mpy(3),h,lunoq(13),ldosq(13),jaq(13) 

39     format(4/,A80)
40     format(6X,6(X,I5))
		  													 
!Inizialización previa de variables:
!-----------------------------------
lunoq=(/0.0651301029002_8, 0.8697397941956_8, 0.0651301029002_8, 0.3128654960049_8, 0.6384441885698_8, 0.0486903154253_8, 	 &
0.6384441885698_8, 0.3128654960049_8, 0.0486903154253_8, 0.2603459660790_8, 0.4793080678419_8, 0.2603459660790_8, 1.0_8/3.0_8/)		   
ldosq=(/0.0651301029002_8, 0.0651301029002_8, 0.8697397941956_8, 0.0486903154253_8, 0.3128654960049_8, 0.6384441885698_8, 	 &
0.0486903154253_8, 0.6384441885698_8, 0.3128654960049_8, 0.2603459660790_8, 0.2603459660790_8, 0.4793080678419_8, 1.0_8/3.0_8/)
jaq=(/0.0266736178004_8, 0.0266736178004_8, 0.0266736178004_8, 0.03855688044515_8, 0.03855688044515_8, 0.03855688044515_8, 	 &
0.03855688044515_8, 0.03855688044515_8, 0.03855688044515_8, 0.0878076287166_8, 0.0878076287166_8, 0.0878076287166_8, 		 &
-0.0747850222335_8/)
lsic=1.0_8
aa=1.0_8    	
	  
!Cálculo de los valores para cada elemento
!-----------------------------------------
open(unit=1,file='C:\malla.txt',status='old')
read(1,39)ac
 do u=1,j
 read(1,40) n(1),n(4),n(2),n(5),n(3),n(6)
 
	 xa=x(n(1))-x(n(3)) 
	 xb=x(n(2))-x(n(3)) 
	 ya=y(n(1))-y(n(3))
	 yb=y(n(2))-y(n(3))  
	 jac=abs(xa*yb-xb*ya)	

     !Cálculo del diametro de la circunferencia inscrita
	 a=sqrt((x(n(1))-x(n(2)))**2.0_8+(y(n(1))-y(n(2)))**2.0_8)
     b=sqrt((x(n(2))-x(n(3)))**2.0_8+(y(n(2))-y(n(3)))**2.0_8)
     c=sqrt((x(n(3))-x(n(1)))**2.0_8+(y(n(3))-y(n(1)))**2.0_8)
     h=2.0_8*jac/(a+b+c)
  
	 !Cálculo de los parámetros de estabilización para elementos que cumplen la condición LBB. 
	 !Aquí sólo es necesario calcular el parámetro 'ccx'.
     ccx=(h**2.0_8)/lsic

	 !Cálculo del módulo de la velocidad en cada nodo esquina y de la distancia elemental para el cálculo del radio hidráulico que se estima 
	 !como el diámetro del círculo de área igual a la del elemento.
	 a=sqrt(vv(n(1))**2+vv(i+n(1))**2)
	 b=sqrt(vv(n(2))**2+vv(i+n(2))**2)
	 c=sqrt(vv(n(3))**2+vv(i+n(3))**2)
	 dist=sqrt((2.0_8*jac)/(3.14159_8))
	 
	 Mui=-(vv(n(1))+vv(n(2))+vv(n(3)))/9.0_8+(vv(n(4))+vv(n(5))+vv(n(6)))*4.0_8/9.0_8   
	 Mvi=-(vv(n(1)+i)+vv(n(2)+i)+vv(n(3)+i))/9.0_8+(vv(n(4)+i)+vv(n(5)+i)+vv(n(6)+i))*4.0_8/9.0_8
	 Mhi=(vv(n(1)+2*i)+vv(n(2)+2*i)+vv(n(3)+2*i))/3.0_8
	 maning=(ma(n(1))+ma(n(2))+ma(n(3)))/3.0_8
	 
	 if (((a.eq.0.0).and.(b.eq.0.0)).or.((a.eq.0.0).and.(c.eq.0.0)).or.((b.eq.0.0).and.(c.eq.0.0))) then
	  !Elementos pegados a un contorno con velocidades nulas. Se consideran sólo elementos con dos nodos apoyados en la pared.
	  s=((maning**2.0_8))/(((dist*Mhi)/(dist+Mhi))**(4.0_8/3.0_8))	  
	 else 
	  !Se consideran elementos con un nodo apoyado en la pared o con ningún nodo apoyado.
	  s=((maning**2.0_8))/(Mhi**(4.0_8/3.0_8)) 
	 endif
	 
	 Mpx(1)=(y(n(2))-y(n(3)))/jac
	 Mpx(2)=(y(n(3))-y(n(1)))/jac
	 Mpx(3)=(y(n(1))-y(n(2)))/jac
	 
	 Mpy(1)=(x(n(3))-x(n(2)))/jac
	 Mpy(2)=(x(n(1))-x(n(3)))/jac
	 Mpy(3)=(x(n(2))-x(n(1)))/jac
	  
	 do uu=1,13	      
	 
	 Mi(1)=lunoq(uu)*(2.0_8*lunoq(uu)-1.0_8)
	 Mi(2)=ldosq(uu)*(2.0_8*ldosq(uu)-1.0_8)
	 Mi(3)=(1.0_8-lunoq(uu)-ldosq(uu))*(1.0_8-2.0_8*lunoq(uu)-2.0_8*ldosq(uu))
	 Mi(4)=4.0_8*lunoq(uu)*ldosq(uu)
	 Mi(5)=4.0_8*ldosq(uu)*(1.0_8-lunoq(uu)-ldosq(uu))
	 Mi(6)=4.0_8*lunoq(uu)*(1.0_8-lunoq(uu)-ldosq(uu))

	 Mix(1)=((4.0_8*lunoq(uu)-1.0_8)*yb)/jac
	 Mix(2)=((1.0_8-4.0_8*ldosq(uu))*ya)/jac
	 Mix(3)=((4.0_8*lunoq(uu)+4.0_8*ldosq(uu)-3.0_8)*(yb-ya))/jac
	 Mix(4)=(4.0_8*ldosq(uu)*yb-4.0_8*lunoq(uu)*ya)/jac
	 Mix(5)=((8.0_8*ldosq(uu)+4.0_8*lunoq(uu)-4.0_8)*ya-4.0_8*ldosq(uu)*yb)/jac
	 Mix(6)=((4.0_8-4.0_8*ldosq(uu)-8.0_8*lunoq(uu))*yb+4.0_8*lunoq(uu)*ya)/jac
	
	 Miy(1)=((1.0_8-4.0_8*lunoq(uu))*xb)/jac
	 Miy(2)=((4.0_8*ldosq(uu)-1.0_8)*xa)/jac
	 Miy(3)=((4.0_8*lunoq(uu)+4.0_8*ldosq(uu)-3.0_8)*(xa-xb))/jac
	 Miy(4)=(4.0_8*lunoq(uu)*xa-4.0_8*ldosq(uu)*xb)/jac
	 Miy(5)=(4.0_8*ldosq(uu)*xb+(4.0_8-4.0_8*lunoq(uu)-8.0_8*ldosq(uu))*xa)/jac
	 Miy(6)=((8.0_8*lunoq(uu)+4.0_8*ldosq(uu)-4.0_8)*xb-4.0_8*lunoq(uu)*xa)/jac

	 Mu=dot_product(vv(n),Mi)   
	 Mv=dot_product(vv(n+i),Mi)
	  
	  do uuu=1,6
	  Nx(n(uuu))=Nx(n(uuu))+ccx*((Mu*Mix(uuu)+Mv*Miy(uuu))*Mu*sqrt(Mui**2+Mvi**2)*s)*jaq(uu)*jac 	 
	  Ny(n(uuu))=Ny(n(uuu))+ccx*((Mu*Mix(uuu)+Mv*Miy(uuu))*Mv*sqrt(Mui**2+Mvi**2)*s)*jaq(uu)*jac 
	  enddo
	  do uuu=1,3
	  Nz(n(uuu))=Nz(n(uuu))+aa*ccx*(Mpx(uuu)*Mu*sqrt(Mui**2+Mvi**2)*s+Mpy(uuu)*Mv*sqrt(Mui**2+Mvi**2)*s)*jaq(uu)*jac
	  enddo	  	  	
	 enddo
 enddo
 close(1)
 end

!--------------------------------------------------------------------------------------------------------------------------------------
!Cálculo de otras variables preproceso o postproceso, tras resolver las ecuaciones diferenciales. Se evalúan las funciones en los nodos
!en vez de en los puntos de integración.
!--------------------------------------------------------------------------------------------------------------------------------------
!Subrutina VELOCIDADESSUBTERRANEAS
!Cálculo postproceso de velocidades	subterráneas o velocidades de Darcy
!--------------------------------------------------------------------------------------------------------------------------------------
subroutine velocidadessubterraneas (i,j,x,y,kix,kiy,ag,vv,velx,vely)
use elemental
integer*4, dimension(:),allocatable::s
integer*4 i,j
real*8, dimension(:),allocatable::kxx,kyy,kxy
real*8 x(i),y(i),kix(i),kiy(i),vv(i),velx(i),vely(i),Mpx(3),Mpy(3),ag(i),Mpxx,Mpyy		 

allocate(s(i),kxx(i),kyy(i),kxy(i))

 39     format(4/,A80)
 40     format(6X,6(X,I5))

!Inizialización previa de variables:
!-----------------------------------
do u=1,i
velx(u)=0.0_8
vely(u)=0.0_8
!Con la variable 's' se hará la media de todos los valores calculados en cada nodo.
s(u)=0
enddo
do u=1,i
kxx(u)=kix(u)*(cos(ag(u))**2)+kiy(u)*(sin(ag(u))**2)
kyy(u)=kix(u)*(sin(ag(u))**2)+kiy(u)*(cos(ag(u))**2)
kxy(u)=(kix(u)-kiy(u))*sin(ag(u))*cos(ag(u))
enddo

!Cálculo de los valores para cada nodo de cada elemento
!------------------------------------------------------
open(unit=3,file='C:\mallasub.txt',status='old')
read(3,39)ac	 
 do u=1,j
     read(3,40) n(1),n(4),n(2),n(5),n(3),n(6)
  
	 xa=x(n(1))-x(n(3)) 
	 xb=x(n(2))-x(n(3)) 
	 ya=y(n(1))-y(n(3))
	 yb=y(n(2))-y(n(3))  	 
	 jac=abs(xa*yb-xb*ya)
   	
	 Mpx(1)=(y(n(2))-y(n(3)))/jac
	 Mpx(2)=(y(n(3))-y(n(1)))/jac
	 Mpx(3)=(y(n(1))-y(n(2)))/jac
	 
	 Mpy(1)=(x(n(3))-x(n(2)))/jac
	 Mpy(2)=(x(n(1))-x(n(3)))/jac
	 Mpy(3)=(x(n(2))-x(n(1)))/jac
	 
	 Mpxx=Mpx(1)*vv(n(1))+Mpx(2)*vv(n(2))+Mpx(3)*vv(n(3))
	 Mpyy=Mpy(1)*vv(n(1))+Mpy(2)*vv(n(2))+Mpy(3)*vv(n(3))

	 !En este caso el valor de la función (la derivada del nivel freático) es constante en el elemento dado que el nivel freático es
	 !un plano que pasa por los tres nodos del elemento y su pendiente es constante. Tendrá el mismo valor en todos los puntos.	
	 !Dado que la derivada no es contínua y toma distintos valores en un nodo entre distintos elementos habrá que hacer la media.	 
	 do uu=1,3
	 !Se suman valores sobre los que pueda haber y se lleva la cuenta en 's' el número de coeficientes considerados por nodo.
	 velx(n(uu))=velx(n(uu))+kxx(n(uu))*Mpxx+kxy(n(uu))*Mpyy 
	 vely(n(uu))=vely(n(uu))+kxy(n(uu))*Mpxx+kyy(n(uu))*Mpyy 
	 s(n(uu))=s(n(uu))+1
	 enddo	 	 
 enddo

 !Cálculo de los valores para cada nodo
 !-------------------------------------
 do u=1,i
  !Se hace la media si se ha considerado más de un valor en ese nodo. En otro caso el valor ya está bien calculado.
  if (s(u).ge.1) then
  velx(u)=-velx(u)/s(u)
  vely(u)=-vely(u)/s(u)
  endif
 enddo 
 close(3)
 deallocate(s,kxx,kyy,kxy)
end
	 
!-------------------------------------------------------------------------------------------------------------------------------
!Subrutina PENDIENTES
!Cálculo de pendientes medias del terreno (en direcciones x e y) multiplicadas por una constante 'mod' que define el usuario
!fuera de la subrutina.   !cñ
!-------------------------------------------------------------------------------------------------------------------------------
subroutine pendientes (i,j,x,y,z,mod,vib)
use elemental
integer*4, dimension(:),allocatable::s
integer*4 i,j
real*8 x(i),y(i),z(i),vib(2*i),Mix(6),Miy(6),Mpxx,Mpyy,mod 		 

allocate(s(i))

 39     format(4/,A80)
 40     format(6X,6(X,I5))

!Inizialización previa de variables:
!-----------------------------------
!En 'lunot,ldost', que son variables globales, van los nodos del elemento (coordenadas naturales).
!Ciertos valores son nulos dado que su dimensión es 7 y sólo se tienen 6 nodos.
lunot=(/1.0_8, 0.0_8, 0.0_8, 0.5_8, 0.0_8, 0.5_8, 0.0_8/)		   
ldost=(/0.0_8, 1.0_8, 0.0_8, 0.5_8, 0.5_8, 0.0_8, 0.0_8/)
do u=1,i
vib(u)=0.0_8
vib(i+u)=0.0_8
s(u)=0
enddo

!Cálculo de los valores para cada nodo de cada elemento
!------------------------------------------------------
open(unit=1,file='C:\malla.txt',status='old')
read(1,39)ac	 
 do u=1,j
     read(1,40) n(1),n(4),n(2),n(5),n(3),n(6)
  
	 xa=x(n(1))-x(n(3)) 
	 xb=x(n(2))-x(n(3)) 
	 ya=y(n(1))-y(n(3))
	 yb=y(n(2))-y(n(3))  	 
	 jac=abs(xa*yb-xb*ya)
   	
	 do uu=1,6
	 !En este caso el valor de la función (la derivada de la cota del terreno) no es constante en el elemento al evaluarse la 
	 !cota del terreno con funciones cuadráticas. Por tanto se evalúa en cada nodo y es necesario tener 'lunot,ldost'.  
	 Mix(1)=((4.0_8*lunot(uu)-1.0_8)*yb)/jac
	 Mix(2)=((1.0_8-4.0_8*ldost(uu))*ya)/jac
	 Mix(3)=((4.0_8*lunot(uu)+4.0_8*ldost(uu)-3.0_8)*(yb-ya))/jac
	 Mix(4)=(4.0_8*ldost(uu)*yb-4.0_8*lunot(uu)*ya)/jac
	 Mix(5)=((8.0_8*ldost(uu)+4.0_8*lunot(uu)-4.0_8)*ya-4.0_8*ldost(uu)*yb)/jac
	 Mix(6)=((4.0_8-4.0_8*ldost(uu)-8.0_8*lunot(uu))*yb+4.0_8*lunot(uu)*ya)/jac
	
	 Miy(1)=((1.0_8-4.0_8*lunot(uu))*xb)/jac
	 Miy(2)=((4.0_8*ldost(uu)-1.0_8)*xa)/jac
	 Miy(3)=((4.0_8*lunot(uu)+4.0_8*ldost(uu)-3.0_8)*(xa-xb))/jac
	 Miy(4)=(4.0_8*lunot(uu)*xa-4.0_8*ldost(uu)*xb)/jac
	 Miy(5)=(4.0_8*ldost(uu)*xb+(4.0_8-4.0_8*lunot(uu)-8.0_8*ldost(uu))*xa)/jac
	 Miy(6)=((8.0_8*lunot(uu)+4.0_8*ldost(uu)-4.0_8)*xb-4.0_8*lunot(uu)*xa)/jac
	 
	 Mpxx=dot_product(Mix,z(n))
	 Mpyy=dot_product(Miy,z(n))  
	 
	 !Dado que la derivada sólo es contínua en la dirección de los lados del elemento habrá que hacer la media.
	 vib(n(uu))=vib(n(uu))+mod*Mpxx 
	 vib(i+n(uu))=vib(i+n(uu))+mod*Mpyy 
	 s(n(uu))=s(n(uu))+1
	 enddo	 	 
 enddo

 !Cálculo de los valores para cada nodo
 !-------------------------------------
 do u=1,i
 if (s(u).ge.1) then
 vib(u)=-vib(u)/s(u)
 vib(i+u)=-vib(i+u)/s(u)
 endif
 enddo 
 close(1)
 deallocate(s)
end

!---------------------------------------------------------------------------------------------------------------------------
!Subrutina VORTICIDAD
!Cálculo postproceso de la vorticidad (en 1/s), de las tensiones en ambas direcciones del espacio y de la tensión 
!tangencial(en Pa=N/m2=kg/(m*s2)).
!---------------------------------------------------------------------------------------------------------------------------
subroutine vorticidad (i,j,x,y,nu,vv,ten,vor)
use elemental
integer*4, dimension(:),allocatable::s
integer*4 i,j
real*8 x(i),y(i),vv(3*i),ten(3*i),vor(i),Mix(6),Miy(6),Mpxx,Mpyy,Mpxy,Mpyx,Mpp,nu,Mp(3)		 

allocate(s(i))

 39     format(4/,A80)
 40     format(6X,6(X,I5))

!Inizialización previa de variables:
!-----------------------------------
lunot=(/1.0_8, 0.0_8, 0.0_8, 0.5_8, 0.0_8, 0.5_8, 0.0_8/)		   
ldost=(/0.0_8, 1.0_8, 0.0_8, 0.5_8, 0.5_8, 0.0_8, 0.0_8/)
do u=1,i
ten(u)=0.0_8
ten(i+u)=0.0_8
ten(2*i+u)=0.0_8
vor(u)=0.0_8
s(u)=0
enddo

!Cálculo de los valores para cada nodo de cada elemento
!------------------------------------------------------
open(unit=1,file='C:\malla.txt',status='old')
read(1,39)ac	 
 do u=1,j
     read(1,40) n(1),n(4),n(2),n(5),n(3),n(6)
  
     xa=x(n(1))-x(n(3)) 
     xb=x(n(2))-x(n(3)) 
     ya=y(n(1))-y(n(3))
     yb=y(n(2))-y(n(3))  	 
     jac=abs(xa*yb-xb*ya)
   	
     do uu=1,6

     Mp(1)=lunot(uu)
     Mp(2)=ldost(uu)
     Mp(3)=1.0_8-lunot(uu)-ldost(uu)

     Mix(1)=((4.0_8*lunot(uu)-1.0_8)*yb)/jac
     Mix(2)=((1.0_8-4.0_8*ldost(uu))*ya)/jac
     Mix(3)=((4.0_8*lunot(uu)+4.0_8*ldost(uu)-3.0_8)*(yb-ya))/jac
     Mix(4)=(4.0_8*ldost(uu)*yb-4.0_8*lunot(uu)*ya)/jac
     Mix(5)=((8.0_8*ldost(uu)+4.0_8*lunot(uu)-4.0_8)*ya-4.0_8*ldost(uu)*yb)/jac
     Mix(6)=((4.0_8-4.0_8*ldost(uu)-8.0_8*lunot(uu))*yb+4.0_8*lunot(uu)*ya)/jac
	
     Miy(1)=((1.0_8-4.0_8*lunot(uu))*xb)/jac
     Miy(2)=((4.0_8*ldost(uu)-1.0_8)*xa)/jac
     Miy(3)=((4.0_8*lunot(uu)+4.0_8*ldost(uu)-3.0_8)*(xa-xb))/jac
     Miy(4)=(4.0_8*lunot(uu)*xa-4.0_8*ldost(uu)*xb)/jac
     Miy(5)=(4.0_8*ldost(uu)*xb+(4.0_8-4.0_8*lunot(uu)-8.0_8*ldost(uu))*xa)/jac
     Miy(6)=((8.0_8*lunot(uu)+4.0_8*ldost(uu)-4.0_8)*xb-4.0_8*lunot(uu)*xa)/jac
	 
	 !Cálculo de las derivadas de la velocidad en cada nodo.
     Mpxx=dot_product(Mix,vv(n))
     Mpyy=dot_product(Miy,vv(i+n))
     Mpxy=dot_product(Miy,vv(n))
     Mpyx=dot_product(Mix,vv(i+n))  
	 !Cálculo del calado en cada nodo. Sólo habrá un valor en cada nodo (no es una derivada), pero la media se hará de todos modos 
	 !dado que no se calculará de forma independiente el término donde el calado aparece.
	 Mpp=Mp(1)*vv(2*i+n(1))+Mp(2)*vv(2*i+n(2))+Mp(3)*vv(2*i+n(3))
	 
	 !Cálculo de tensiones en direcciones x,y y tensión tangencial. 
	 ten(n(uu))=ten(n(uu))+1000.0_8*(-9.81_8*Mpp+2.0_8*Mpxx*nu) 
	 ten(i+n(uu))=ten(i+n(uu))+1000.0_8*(-9.81_8*Mpp+2.0_8*Mpyy*nu)
	 ten(2*i+n(uu))=ten(2*i+n(uu))+1000.0_8*nu*(Mpxy+Mpyx) 
	 !Cálculo de la vorticidad.
	 vor(n(uu))=vor(n(uu))+Mpyx-Mpxy
	 s(n(uu))=s(n(uu))+1
	 enddo	 	 
 enddo

 !Cálculo de los valores para cada nodo
 !-------------------------------------
 do u=1,i
 if (s(u).gt.1) then
 ten(u)=ten(u)/s(u)
 ten(i+u)=ten(i+u)/s(u)
 ten(2*i+u)=ten(2*i+u)/s(u)
 vor(u)=vor(u)/s(u)
 endif
 enddo 
 close(1)
 deallocate(s)
end

!---------------------------------------------------------------------------------------------------------------------------------------
!Subrutinas para la construcción de la matriz del sistema.
!---------------------------------------------------------------------------------------------------------------------------------------
!Subrutina SUMA
!Se llama a esta subrutina cada vez que se forma una matriz elemental. Con el almacenamiento mediante vectores sólo es necesario 
!dimensionar las matrices elementales de dimensión seis por seis aquí llamadas 'Nn'. Es necesario guardar los valores enteros n(6) que 
!llevan el número de los 6 nodos de cada elemento para poder colocar en su posición los coeficientes de cada matriz elemental. También 
!los valores enteros 'e' y 'ee' que sitúan las cajas en el sistema en su posición. Aquí sa es el vector donde van los coeficientes 
!considerados e 'ita' es el vector puntero que lleva la posición de estos coeficientes (se utiliza formato tipo MSR donde los 
!coeficientes no están ordenados por columnas, donde pueden aparecer coeficientes nulos, coeficientes de la diagonal o coeficientes de 
!la misma posición). 'inc' es el tamaño total del sistema creado inicialmente, antes de usar la subrutina reducción.
!---------------------------------------------------------------------------------------------------------------------------------------
subroutine suma(Nn,n,e,ee,inc,uj,ui,po)
use allocatacion
integer*4 n(6),r,s,rr,ss,e,ee,inc,uj,ui,po(inc)  
real*8 Nn(uj,ui)

!Escritura de coeficientes uno tras otro
!---------------------------------------
!Las primeras 'inc+1' componentes de 'ita' ya han sido referenciadas en la subrutinas dimvectas o dimvectsb.
do r=1,uj
rr=n(r)+e
do s=1,ui
ss=n(s)+ee 
sa(po(rr))=sa(po(rr))+Nn(r,s)
ita(po(rr))=ss
po(rr)=po(rr)+1
enddo
enddo
end

!------------------------------------------------------------------------------------------------------------------------
!Subrutina ORDEN
!Ordenamiento de los coeficientes almacenados y generación del formato MSR para los vectores 'ita,sa'.
!Habrá que tener este formato si se quiere resolver el sistema o si se quiere hacer un producto matriz-vector.
!------------------------------------------------------------------------------------------------------------------------
subroutine orden (inc,k)
use allocatacion
integer*4 k,r,s,inc,u  
real*8 kr

!Ordenamiento de de los coeficientes de cada fila 
!------------------------------------------------
!Con este proceso aparecerán coeficientes con la misma posición uno al lado del otro en el vector.
do u=1,inc
do r=ita(u),ita(u+1)-1	  
do s=r,ita(u+1)-1
 if (ita(s).lt.ita(r)) then
 k=ita(r)
 kr=sa(r)
 ita(r)=ita(s)
 sa(r)=sa(s)
 ita(s)=k
 sa(s)=kr
 endif
enddo
enddo
enddo

!Ensamblamiento y eliminación de los posibles ceros y del término de la diagonal
!-------------------------------------------------------------------------------
!Se trata de procesos que cambian 'ndim'.
!Primero se ensambla sumando coeficientes con la misma posición. Se quitan los posibles ceros, los introducidos inicialmente 
!y los que aparecen tras el ensamblaje. También se quita el término de la diagonal (éste se almacena en su posición).  
k=0
do u=1,inc 
!En primer lugar se escribe el primer coeficiente de cada fila teniendo en cuenta la posición donde está.
if ((ita(u)+k).ne.ita(u+1)) then 
!Si no hay coeficientes en la fila (si 'ita(u)' original es igual a 'ita(u+1)') no se accede a ella
ita(ita(u))=ita(ita(u)+k)
sa(ita(u))=sa(ita(u)+k)
!En segundo lugar se escriben los siguientes coeficientes de la fila ensamblándolos cuando sea necesario.
do r=ita(u)+k,ita(u+1)-2
 if (ita(r+1).eq.ita(r)) then
  !Se debe ensamblar el término 'r+1'.
  !Se espera a evaluar el coeficiente con el siguiente paso del bucle si se ensambla. Sería posible ensamblar de nuevo.	 
  sa(r-k)=sa(r-k)+sa(r+1)
  k=k+1
 else
  !No se ensambla el término 'r+1'.
  !En tercer lugar se evalúan todos los coeficientes escritos para su posible eliminación.
  !Se evalúa el término 'r-k' donde se habrán ensamblado coeficientes en caso de haber sido necesario. 
  if (ita(r-k).eq.u) then
  !Se evalúa si es un coeficiente de la diagonal.
  sa(u)=sa(r-k) 
  k=k+1
  elseif (sa(r-k).eq.0.0_8) then
  !Se evalúa si es un coeficiente nulo.
  k=k+1
  endif 
 !Se escribe el coeficiente 'r+1'. Éste no se evaluará dentro de este bucle si se trata del último coeficiente de la fila.
 sa(r+1-k)=sa(r+1)
 ita(r+1-k)=ita(r+1)
 endif
enddo
!Análisis del último término escrito en esa fila. Procedimiento válido si sólo hay un coeficiente en la fila, en cuyo caso
!no se ha entrado en el bucle anterior pero se ha escrito inicialmente el primer y único coeficiente.
r=ita(u+1)-1
 if (ita(r-k).eq.u) then
 !Se evalúa si es un coeficiente de la diagonal.
 sa(u)=sa(r-k) 
 k=k+1
 elseif (sa(r-k).eq.0.0_8) then
 !Se evalúa si es un coeficiente nulo.
 k=k+1
 endif
endif
!Modificación de las primeras 'inc+1' componentes del vector 'ita'.
ita(u+1)=ita(u+1)-k
enddo
end

!------------------------------------------------------------------------------------------------------------------------------------------
!Imposición de CC reduciendo el orden del sistema.
!------------------------------------------------------------------------------------------------------------------------------------------
!Subrutina REDUCCIONDELSISTEMA
!Se reducirá el sistema para flujo superficial (a excepción del que aparece con el método de Newton) y para flujo subterráneo al imponer 
!las condiciones de contorno. 
!Las condiciones de contorno vienen desde el programa principal guardadas en 'vic' si no existe un contorno móvil. Éstas habrán sido leídas 
!del fichero 'malla.txt', del fichero 'mallasub.txt'. En otro caso las condiciones vienen desde las subrutinas que gestionan el movimiento 
!de los contornos móviles. Éstas habrán sido leídas del fichero 'malla.txt', del fichero 'mallasub.txt' o bien dadas en ejecución.
!Los contadores 'i,j' son el número de nodos y de elementos que vienen desde la subrutina aguassomeras o la subrutina aguassubterranea 
!según se resuelva el sistema para flujo superficial o subterráneo. Aquí 'inc' será diferente para cada sistema.
!Se modificará la matriz almacenada en 'isa,sa' una vez que se le de formato MSR (y sea ensamblada) y el término independiente almacenado 
!en 'vector' que ya está ensamblado. Se llamará a la subrutina reduccion para modificar la matriz manteniendo el formato MSR.  
!-----------------------------------------------------------------------------------------------------------------------------------------
subroutine reducciondelsistema(i,c,vector,vic,v,vb,inc,bcg,ndim)
use allocatacion
integer*4, dimension(:),allocatable::ck
integer*4 i,u,uu,v(i),c,ndim,inc,w,k
real*8, dimension(:),allocatable::vi
real*8 vector(inc),vic(inc),vb(inc)
character bcg*2	   

allocate(ck(3*i),vi(inc))

!Escritura de la matriz del sistema en formato MSR (Modified Sparse Row)
!-----------------------------------------------------------------------
!Generación de vectores con formato MSR.
!Se tiene que 'inc=3*i' para las ecuaciones de aguas someras o de Navier-Stokes 2D e 'inc=i' para la ecuación de agua subterránea.
call orden (inc,k)
ndim=ndim-k

!Casos en que se reduce la matriz y el término independiente:
!Primer caso de reducción: Los pesos de la ecuación de continuidad (en ambos modelos) sólo están discretizados para los nodos esquina y 
!el calado en las dinámicas (modelo superficial) o el nivel freático en la ecuación de continuidad (modelo subterráneo) sólo están discretizados 
!en nodos esquina. Por ello es necesario eliminar las filas y columnas de las cajas D, E, B, Asx, Asy, Ns, N,... (sólo habrá coeficientes nulos) 
!y los respectivos coeficientes del término independiente. Se reducirá la dimensión del sistema de '3*i' a '2*i+nodos esquina' en el modelo 
!superficial y de 'i' a 'nodos esquina' en el subterráneo. 
!Dar CC nulas genera el sistema que se debería usar con la discretización menor (elimina fila y columna simplemente).

!Segundo caso de reducción: Se produce al imponer las condiciones de contorno (CC) guardadas en el vector 'vic', pasando el sistema de tener 
!la dimensión señalada a tener '2*i+(nodos esquina)-CC' en el modelo superficial y '(nodos esquina)-CC' en el modelo subterránea, donde la 
!referencia de los nodos esquina está en el vector v(inc). Si no se usa toda la malla habrá más condiciones que las definidas en el contorno y 
!contorno móvil. Esto es, la parte de la malla que no se usa con las ecuaciones superficiales (con nodos secos para el modelo conjunto o 
!superficial) o la parte de la malla que no se usa con las ecuaciones subterráneas (modelo conjunto) no se consideran y para ello se han dado 
!CC=0 en estos nodos. 
!Las CC nulas permiten generar el sistema que se debería usar al no consideran las partes de la malla que tienen esa condición.

!Imposición de condiciones de contorno sobre el sistema
!------------------------------------------------------
!Se modifica el vector 'vic' de condiciones de contorno añadiendo las CC nulas para el primer caso y así se efectúa la reducción conjuntamente.
c=0
w=0
do u=1,i
 if (v(u).eq.1) then
 vic(u+inc-i)=0.0_8
 endif
enddo 

!Se modifica el término independiente con las CC almacenadas en 'vic' y se reduce el orden del término independiente (se eliminan filas).
do u=1,inc 
 if (vic(u).eq.sqrt(2.0_8)) then
 vi(u)=0.0_8
 c=c+1
 else
 vi(u)=vic(u)
 endif
enddo

do u=1,c
do while (vic(u+w).ne.sqrt(2.0_8)) 
 w=w+1
enddo
 vector(u)=vector(u+w)
 vb(u)=vb(u+w)
 do uu=ita(u+w),ita(u+w+1)-1
 vector(u)=vector(u)-sa(uu)*vi(ita(uu))
 enddo 
enddo

!Se pasa de una matriz cuadrada de dimensión inc a otra cuadrada con la dimensión señalada en cada caso, trabajando sobre los vectores 
!'ita,sa' en que está escrita la matriz del sistema.
call reduccion(ndim,vic,inc)
	  
!El sistema a resolver ahora es de dimension 'c' y la dimensión de 'ita,sa' es 'ndim' que se puede calcular como 'ita(c+1)-1' (formato MSR)
write(6,*)' '
write(6,*)'Longitud del vector matriz dispersa:',ndim

!Escritura de la matriz del sistema en formato CSC (Compressed Sparse Column)
!----------------------------------------------------------------------------
!La dimensión de los vectores que almacenan la matriz será 'ndim-1' para dos de los vectores que se puede calcular como 'ica(c+1)-1'  
!(formato CSC) y será 'c+1' para el tercero.			  
if (bcg.eq.'no') then
 !Se deja el formato MSR si se va a calcular con precondicionador diagonal. Los vectores 'cia,ca' no son necesarios.
 deallocate(cia,ca)
else
 !Se pasa a formato CSC si se va a calcular con precondicionador LU.
 allocate(cja(c+1))
 !Inicialización de variables.
 do u=1,3*i
 ck(u)=0
 enddo
 cja(1)=1
 do u=1,c
 cja(u+1)=0
 enddo
 !Cálculo del número de coeficientes por columna. 
 do u=c+2,ndim
 cja(ita(u)+1)=cja(ita(u)+1)+1
 enddo
 !Configuración final de 'cja' (considerando el coeficiente de la diagonal) y escritura de la diagonal en 'cia,ca'.
 !Además se tiene el contador 'ck'.
 do u=1,c						
 cja(u+1)=cja(u+1)+cja(u)+1
 ca(cja(u))=sa(u)
 cia(cja(u))=u
 ck(u)=cja(u)+1
 enddo
 !Escritura del resto de coeficientes en 'cia,ca'.
 do u=1,c										   
  do w=ita(u),ita(u+1)-1
  !En 'w' está la columna y en 'ck(w)' la posición preparada para ese coeficiente.
  cia(ck(ita(w)))=u
  ca(ck(ita(w)))=sa(w)
  ck(ita(w))=ck(ita(w))+1
  enddo
 enddo
 deallocate(ita,sa)
endif
!Ahora la dimensión de 'ca,cia' es 'ndim-1' que se puede calcular como 'cja(c+1)-1'.
!Los vectores 'ita,sa' no son necesarios.
deallocate(ck,vi) 
end

!-----------------------------------------------------------------------------------------------------------------------------------------------
!Subrutina REDUCCIÓN
!En esta subrutina se eliminan todas las filas y columnas referenciadas en 'vic' a la vez. Si existe valor en 'vic(u)' se eliminará
!la fila y la columna 'u' del sistema existente en el momento de llamar a esta subrutina. 
!Se operará directamente sobre los vectores 'sa,ita' que son respectivamente el vector donde van los coeficientes no nulos y el vector puntero 
!(formato MSR). 'ndim' es la dimensión de los vectores 'ita,sa' y 'iya,ya' guardan las primeras 'inc+1' componentes originales de 'ita,sa'.
!'lo' es el número de filas-columnas a eliminar.
!'vac' es un vector que lleva en cada posición 'u' el número de filas-columnas a eliminar con valor menor o igual que 'u' (considera 
!también la eliminación de la fila-columna 'u' en ese número). 
!'wo' llevará el número de coeficientes no nulos existentes en cada una de las filas a eliminar, sin considerar el coeficiente de la diagonal. 
!Incluye coeficientes correspondientes a otras columnas a eliminar.	Si la fila 'u' se debe eliminar, se deberán eliminar 'wo(u)' coeficientes. 
!----------------------------------------------------------------------------------------------------------------------------------------------
subroutine reduccion(ndim,vic,inc)    
use allocatacion
integer*4, dimension(:),allocatable::vac,wo,iya
integer*4 u,v,w,ndim,inc,lo,le,li 
real*8, dimension(:),allocatable::ya 
real*8 vic(inc)

allocate(vac(inc),wo(inc),iya(inc+1),ya(inc+1))

!Inizialización previa de variables
!----------------------------------
w=0
lo=0
do u=1,inc
 if (vic(u).ne.sqrt(2.0_8)) then
 lo=lo+1
 wo(u)=ita(u+1)-ita(u)
 else
 wo(u)=0
 endif
vac(u)=lo
enddo
do u=1,inc+1
iya(u)=ita(u)
ya(u)=sa(u)
enddo

!Se eliminan las filas 'u' del sistema si 'vic(u)' tiene valor (si es sqrt(2.0), no tiene valor)
!-----------------------------------------------------------------------------------------------
!Se sobreescriben coeficientes durante este proceso.
le=lo
do u=1,inc-lo		
 !'u' serán las filas a considerar (las que no serán eliminadas) y 'w' las que se consideran.
 do while (vic(u+w).ne.sqrt(2.0_8))
 !En 'le' (variable en el bucle) se considera la eliminación los 'wo' coeficientes que se tienen hasta la fila 'u+w' (variable) y los 
 !coeficientes a eliminar en las primeras 'inc+1' componentes de 'ita,sa' originales (constante='lo').
 le=le+wo(w+u)
 w=w+1
 enddo
 !Se colocan los coeficientes (no nulos) de las filas a considerar desde la posición 'inc+2-lo'. Es necesario que 'le' acumule el número 
 !de coeficientes almacenados en 'wo' que ya se han eliminado ya que 'le' será la distancia con que serán traspasados los coeficientes. 
 !Se escriben todos los coeficientes de las filas a considerar de forma contínua. De momento tendrán todos sus coeficientes, incluyendo 
 !aquéllos de las columnas a eliminar.  
 do v=iya(u+w)-le,iya(u+w+1)-1-le
 ita(v)=ita(v+le)
 sa(v)=sa(v+le)
 enddo
!Se va a referenciar la nueva situación de los coeficientes modificando las primeras 'inc+1-lo' componentes del vector 'ita'. En ellas se 
!referenciará la posición del primer coefiente (no nulo) de todas las filas a considerar. Igualmente, se tendrán en cuenta todos sus 
!coeficientes, incluyendo aquéllos de las columnas a eliminar. Así, 'ita(u)' siempre llevará la posición inicial de la fila 'u' considerada. 
ita(u)=iya(u+w)-le
sa(u)=ya(u+w)
enddo
!Para u=inc+1-lo, 'u+w' solo es mayor que 'inc' si la última fila-columna del sistema no se elimina ('w=lo'). En este caso no hará falta 
!recalcular 'le' ni 'w' y no se entrará en el bucle. En otro caso se eliminará la última o últimas filas-columnas y se entrará en el bucle.
u=inc+1-lo
do while ((u+w).le.inc)	 
le=le+wo(w+u)
w=w+1
enddo
ita(u)=iya(u+w)-le
sa(u)=0.0_8  

!Se eliminan las columnas 'u' del sistema si 'vic(u)' tiene valor 
!----------------------------------------------------------------
!'lo' continene el número de filas-columnas a eliminar o coeficientes eliminados en las primeras 'inc+1' posiciones de 'ita,sa' original.
!'le' contiene el número de coeficientes eliminados hasta ahora ('le-lo' son los coeficientes eliminados en las posiciones desde 'inc+2' a 
!'ndim' de 'ita,sa' originales).
!Ahora 'ita,sa' referencian a un sistema rectangular incx(inc-lo) con 'inc' columnas dado que no se han eliminado los coeficientes 
!correspondientes a las columnas a eliminar en las filas que se consideran.
!A continuación se eliminan todas las columnas ('li' coeficientes) para formar un sistema (inc-lo)x(inc-lo) de forma que sólo se pase una vez 
!por cada coeficiente de los vectores 'ita,sa'). Desde el primer coeficiente que se elimina, todos los coeficientes siguientes deben ser 
!reposicionados (a una distancia). De nuevo se sobreescriben coeficientes durante este segundo proceso.
li=0
do u=1,inc-lo										   
w=ita(u)	
 do while ((w+li).ne.ita(u+1))	  
  !Se evalúa una fila.
  if (vic(ita(w+li)).ne.sqrt(2.0_8)) then
  !Se evalúa si el coeficiente pertenece a cualquiera de las columnas a eliminar y en tal caso se modifica 'li'.
  li=li+1
  else
  !En cada fila se reposicionan los coeficientes y la referencia de su posición con 'vac' (sobreescribiendo sobre los coeficientes previos 
  !si pertenecen a columnas a eliminar).
  sa(w)=sa(w+li)	 	   
  ita(w)=ita(w+li)-vac(ita(w+li))		 	    
  w=w+1
  endif
 enddo
!Modificación de las primeras 'inc+1-lo' componentes del vector 'ita'.
ita(u+1)=ita(u+1)-li
enddo

!Cálculo de la dimensión final de los vectores 'ita,sa'
!------------------------------------------------------
!'li' contiene el número de coeficientes eliminados durante el segundo proceso. En total se habrán eliminado 'le+li' coeficientes de los 
!vectores 'ita,sa' originales.
ndim=ndim-le-li
deallocate(vac,wo,iya,ya)
end

!----------------------------------------------------------------------------------------------------------------------------------------
!Método PBCG con precondicionamiento diagonal - Obtenida del libro Fortran recipes (chap.2). Autor:  Greenbaum, Anne, (Courant Institute)
!Precondicionador seleccionado, y código reprogramado con diferentes comandos y adaptado para .f95.
!----------------------------------------------------------------------------------------------------------------------------------------
!Subrutina GRADIENTESBICONJUGADOS
!En esta subrutina se aplica método PBCG con precondicionamiento diagonal (también es posible aplicar el método BCG sin 
!precondicionamiento). Se resuelve el sistema lineal Ax=b. Si A es definida positiva y simétrica se aplica el método CG con más 
!operaciones de las necesarias.
!'i' es el número de nodos de la malla, 'c' es la dimension de la matriz del sistema, 'vec' es el vector de terminos independientes del 
!sistema, 'er' indica la posibilidad de que el error esté mal calculado, 'err' es el error estimado, y 'r,s,k' son contadores. 
!Las iteraciones se paran cuando el valor 'abs(Ax-b)/abs(b)' es menor que 'tol' (también es posible utilizar otro test de convergencia
!modificando 'itol'). El número máximo de iteraciones permitido será '15*c' y se toma 'tol=1e-10'. Si 'tol>err' y 'iter>15*c' se 
!obtendrá la solución del sistema en cualquier caso (con mayor error).
!----------------------------------------------------------------------------------------------------------------------------------------
subroutine gradientesbiconjugados(vec,c,x,nonzero)	   			   
use allocatacion
integer*4 epsil
parameter(epsil=1.d-14)
integer*4 c,r,s,iter,itol,er
real*8, dimension(:),allocatable::p,pp,ra,rr,z,zz,di
real*8 err,tol,x(c),vec(c),he,ak,akden,bk,bkden,bknum,bnrm,dxnrm,xnrm,znrm,snrm 
logical nonzero	 
character tiempo*8

allocate(p(c),pp(c),ra(c),rr(c),z(c),zz(c),di(c))

302   format (A7,I8,5X,A7,E11.4E2)

!Se muestra por pantalla la hora (con formato hh:mm:ss) cuando se entra aquí
!---------------------------------------------------------------------------
call time(tiempo)
write(6,*) 'Hora:',tiempo

!Inicialización previa de variables
!----------------------------------
!Se toma 'itol=1'.
!Si 'itol=1' las iteraciones se paran cuando el valor 'abs(Ax-b)/abs(b)' es menor que 'tol'.
!Si 'itol=2' las iteraciones se paran cuando el valor 'abs(inv(P)(Ax-b))/abs(inv(P)b)' es menor que 'tol'. P es el precondicionador tal
!que P= matrizsimilar(A), y será mejor tanto más similar sea.
!Si 'itol=3' la subrutina usa su propio estimador del error en 'x', y las iteraciones paran cuando su magnitud dividida por la magnitud 
!de 'x' es menor que 'tol'. 
!Si 'itol=4' se aplica lo mismo que si 'itol=3' excepto por el hecho de que la mayor componente del error en valor absoluto y la mayor 
!componente de 'x' son utilizadas en vez de la magnitud del vector.  	  
itol=1
tol=1e-10
err=1.0_8
iter=0											   
he=0.0_8
znrm=1.0_8
dxnrm=0.0_8
xnrm=0.0_8

!Precondicionamiento
!-------------------
!Sin precondicionamiento.
!Buenos resultados en todos los casos aunque a veces hay inestabilidades del método. Se necesitan muchas iteraciones. 
!do r=1,c 
!di(r)=1.0_8     
!enddo

!Con precondicionamiento (precondicionador diagonal) o sin precondicionamiento.
!Si no hay coeficientes nulos en la diagonal de la matriz del sistema se usa la matriz diagonal P con la diagonal de dicha matriz (se 
!precondiciona con la inversa). Ocurre en GW, SW no estacionario con Picard, NS, SW si se usa Newton o estabilización.
!En el GW se mantiene o se reduce el nº iteraciones. En el SW no estacionario parece que no va bien. Para el método de Newton no va bien 
!(tal vez porque los coeficientes de la diagonal en la 3ª ec. son muy pequeños). Para estabilización siempre va mejor en caso estacionario 
!que si no se usa precondicionamiento. 
!En otro caso se utilizará la matriz identidad y se resuelve sin precondicionamiento (caso anterior). Ocurre en NS, SW estacionario 
!(si no se usa Newton).
!do r=1,c
!di(r)=sa(r)
!enddo
!if (minval(abs(di)).eq.0.0_8) then    
! do r=1,c 
! di(r)=1.0_8     
! enddo 
!endif

!Con precondicionamiento (precondicionador diagonal). 
!Se usa la matriz diagonal P donde los coeficientes de la diagonal son la norma de los coeficientes de cada fila (se precondiciona con la 
!inversa). Es válido en todos los casos (haya o no ceros en la diagonal) como ocurre con el método PCGBLU.
!Buenos resultados en todos los casos aunque a veces hay inestabilidades del método. Consigue reducir el número de iteraciones.
do r=1,c
di(r)=0.0_8
 do s=ita(r),ita(r+1)-1
 di(r)=di(r)+sa(s)**2
 enddo
 di(r)=sqrt(di(r)+sa(r)**2)
enddo

!Se calcula el residual inicial dependiendo de la aproximación inicial dada
!--------------------------------------------------------------------------
!La aproximación inicial será o bien un vector con sus coeficientes nulos o un vector con la última solución del sistema obtenida. 
  if (nonzero) then
   !Si el valor inicial 'x' no es nulo. 
   call dsprsax(x,ra,c)
   !Se realiza un producto matriz-vector con las subrutinas dsprsax.
   do s=1,c
   ra(s)=vec(s)-ra(s)
   rr(s)=ra(s)
   enddo  
  else
   !Si el valor inicial 'x' es nulo.
   do s=1,c
   ra(s)=vec(s)
   rr(s)=vec(s)
   enddo
  endif
  !Utilizando aquí 'call dsprsax(ra,rr,c)' se calcularía un residual 'rr' diferente y se tendría el algoritmo de mínimo residual 
  !para A simetrica y no definida positiva (el GMRES para matrices no simétricas es de este tipo). Es una versión del algoritmo CG.

 !Inicialización de variables según el test de convergencia aplicado
 !------------------------------------------------------------------
  if (itol.eq.1) then 
   bnrm=snrm(c,vec,itol)
   !Con asolve se utiliza el precondicionador (P).
   call asolve(c,ra,z,di)
  elseif (itol.eq.2) then
   call asolve(c,vec,z,di)
   bnrm=snrm(c,z,itol)
   call asolve(c,ra,z,di)
  elseif ((itol.eq.3).or.(itol.eq.4)) then
   call asolve(c,vec,z,di)
   bnrm=snrm(c,z,itol)
   call asolve(c,ra,z,di)
   znrm=snrm(c,z,itol)
  else
   pause 'itol debe estar entre 1 y 4'
  endif 
   
   !Comienzo de las iteraciones del método
   !--------------------------------------
   do while (err.gt.tol)
    !Escritura por pantalla de iteracion y del error. 
    if((iter.eq.100).or.(iter.eq.1000).or.(iter.eq.c).or.(iter.eq.7*c).or.(iter.eq.15*c)) then   
	write(6,302) ' iter =',iter,'error =',err    
     if (er.eq.1) then
	 write(6,*)'Error posiblemente inexacto'
	 endif
	 !Parada del bucle para iter=itmax=15*c.
     if(iter.eq.15*c) then
     write(6,*)' '
     write(6,*)'Iteracion maxima permitida y convergencia no alcanzada'
     exit    
     endif 
    endif
	er=0
	iter=iter+1
    !Se debe utilizar la matriz traspuesta de la matriz precondicionamiento (PT) que es igual a la matriz de 
	!precondicionamiento (P) al ser una matriz diagonal. Por tanto se usa la subrutina asolve.
    call asolve(c,rr,zz,di)
    bknum=0.0_8
	 bknum=dot_product(z,rr)
      if (iter.eq.1) then 
       do s=1,c
       p(s)=z(s)
       pp(s)=zz(s)
       enddo
      else
       bk=bknum/bkden
       do s=1,c
       p(s)=bk*p(s)+z(s)
       pp(s)=bk*pp(s)+zz(s)
       enddo
      endif
      !A continuación se calculan los residuales 'ra' y 'rr' para cada iteración.
      bkden=bknum
      call dsprsax(p,z,c)
      akden=0.0_8
       akden=dot_product(z,pp)
      ak=bknum/akden
      !Se realiza un producto matriz-vector con la subrutina dsprstx.
	  call dsprstx(pp,zz,c)
	   do s=1,c
       x(s)=x(s)+ak*p(s)
       ra(s)=ra(s)-ak*z(s)
       rr(s)=rr(s)-ak*zz(s)
       enddo
	  !Con asolve se resuelve P*z=ra y se evalúa el criterio de parada.
      call asolve(c,ra,z,di)
       if (itol.eq.1) then     
        err=snrm(c,ra,itol)/bnrm
	   elseif (itol.eq.2) then
		err=snrm(c,z,itol)/bnrm
       elseif ((itol.eq.3).or.(itol.eq.4)) then
	    he=znrm
        znrm=snrm(c,z,itol)
        if (abs(he-znrm).gt.epsil*znrm) then
         dxnrm=abs(ak)*snrm(c,p,itol)
         err=znrm/abs(he-znrm)*dxnrm 
        else
	     !El error calculado podría ser inexacto. 
         err=znrm/bnrm
		 er=1
         cycle
        endif
         xnrm=snrm(c,x,itol)
        if (err.le.0.5_8*xnrm) then
         err=err/xnrm
        else 
	     !El error calculado podría ser inexacto.                        
         err=znrm/bnrm
		 er=1      
         cycle
        endif
       endif 
   enddo
   
   write(6,*) 'Solucion obtenida en iteracion:',iter
   deallocate(p,pp,ra,rr,z,zz,di)            
end

!---------------------------------------------------------------------------------
!Función SNRM
!Cálculo de la norma para el vector 'sx'.  
!---------------------------------------------------------------------------------
function snrm(c,sx,itol)
integer*4 c,itol,r,isamax
real*8 sx(c),snrm

  if (itol.le.3) then
    !Cálculo de la norma del vector.   
    snrm=sqrt(dot_product(sx,sx))
  else
    !Cálculo de la norma de la mayor componente (es directamente esa componente).
	isamax=1
    do r=1,c
     if (abs(sx(r)).gt.abs(sx(isamax))) then
     isamax=r
     endif
    enddo
    snrm=abs(sx(isamax))
  endif
end          
  
!---------------------------------------------------------------------------------------------------------------------------
!Subrutina ASOLVE
!Precondicionamiento de la matriz asociada al sistema con una matriz diagonal. Se resuelve un sistema P*x=vec con solución 
!x=P(-1)*vec, estando P definido en 'di'. 
!'di' puede contener la diag(I), la diag(A) o bien la norma de los términos de cada fila.
!---------------------------------------------------------------------------------------------------------------------------
subroutine asolve (c,vec,x,di)
integer*4 c,r
real*8 x(c),vec(c),di(c)   

  !Multiplico la inversa de la matriz de precondicionamiento por vec (como en la subrutina dsprsax se multiplica una matriz 
  !por x(r)). La inversa de una matriz diagonal es otra matriz diagonal con diagonal 1/(coeficientes de la diagonal). 
  !x=P(-1)*vec=vec/diag(P)=vec/di
  do r=1,c 
  x(r)=vec(r)/di(r)     
  enddo     
end
	  
!-------------------------------------------------------------------------------------------
!Subrutina DSPRSAX
!Esta subrutina multiplica una matriz dispersa de dimension 'c' almacenada por filas en 'sa' 
!con puntero 'ita' por un vector columna 'x' y lo almacena en 'vec'.
!-------------------------------------------------------------------------------------------
subroutine dsprsax(x,vec,c) 
use allocatacion
integer*4 r,k,c
real*8 vec(c),x(c)
  
  !Producto matriz-vector.
  do r=1,c
   vec(r)=sa(r)*x(r)
   do k=ita(r),ita(r+1)-1
    vec(r)=vec(r)+sa(k)*x(ita(k))
   enddo
  enddo
end  

!--------------------------------------------------------------------------------------------
!Subrutina DSPRSTX
!Esta subrutina multiplica la traspuesta de una matriz dispersa de dimension 'c' almacenada 
!por filas 'sa' con puntero 'ita' por un vector columna 'x' y lo almacena en 'vec'.
!--------------------------------------------------------------------------------------------
subroutine dsprstx(x,vec,c)
use allocatacion
integer*4 r,k,c
real*8 vec(c),x(c)
      
  !Producto matriz-vector.
  do r=1,c
   vec(r)=sa(r)*x(r) 
  enddo 
  do r=1,c 
   do k=ita(r),ita(r+1)-1     
    vec(ita(k))=vec(ita(k))+sa(k)*x(r)
   enddo
  enddo
end  

!---------------------------------------------------------------------------------------------------------------------------------------------
!Método PBCGLU con precondicionamiento LU - Obtenida de la librería SLATEC (SLAP). Autores:  Greenbaum, Anne (Courant Institute), 
!Seager, Mark K. (LLNL). Código reprogramado con diferentes comandos y adaptado para .f95. 
!Se definen varias constantes que dependen de la máquina a través de funciones intrínsecas de Fortran (obviamente 1.d0=1.0_8).
!Magnitud positiva más pequeña = tiny(1.0). En ieee standard (intel 8087 and intel 8086 emulator) su valor es 2.22d-308.
!Magnitud más grande = huge(1.0). En ieee standard su valor es 1.79d308.	  
!Espacio relativo más pequeño = epsilon(1.0)/radix(1.0_8)=epsilon(1.0)/2.0. En ieee standard su valor es 1.11d-16.
!epsilon (1.0) sería el espacio relativo más grande.
!---------------------------------------------------------------------------------------------------------------------------------------------
!Subrutina DSLUBC
!Subrutina para resolver un sistema lineal Ax=b utilizando el método de los gradientes biconjugados con descomposición LU incompleta 
!para formar la matriz de precondicionamiento. Los coeficientes de las matrices que se generarán se guardarán en 'iw y 'w' que se generan 
!aquí y que tendrán dimensiones 'leniw' y 'lenw' respectivamente. 'n' será la dimension de la matriz del sistema. 
!Las iteraciones se paran cuando el valor 'abs(Ax-b)/abs(b)' es menor que 'tol' (también es posible utilizar otro test de convergencia
!modificando 'itl'). 'itl' y 'tl' son las mismas variables que 'itol,tol' en la subrutina gradientesbiconjugados del método PBCG y toman los 
!mismos valores. La subrutina devuelve 'ite' y 'er' como la iteración y error estimado alcanzados. 'itmax' es el nº máximo de iteraciones 
!permitidas y será '15*c'(en la subrutina gradientesbiconjugados lo evaluaba dentro del bucle). Si 'ite>itmax' se obtendrá la solución del 
!sistema en cualquier caso (con mayor error). Se obtendrá 'ite=itmax+1' si 'er>tl'.
!La subrutina devuelve 'ie' con el que se podrá saber el tipo de error que se generó en caso de que haya algún problema.
!'ie=0' si todo fue bien, 'ie=1' si hay insuficiente espacio reservado para 'w,iw', 'ie=2' si no se obtiene convergencia para el número máximo 
!de iteraciones, 'ie=3' si hay error en los datos de entrada de 'n' o 'itol', 'ie=4' si la tolerancia es demasiado estricta en cuyo caso fue
!reseteada a 500*epsilon(1.0)/2.0, 'ie=5' si la matriz de precondicionamiento no es definida positiva. Producto escalar (denominador de bk) 
!(r,z)<0. Se evalúa con una tolerancia, 'ie=6' si la matriz A no es definida positiva. Producto escalar (denominador de ak)  (p,A*p)<0. 
!Se evalúa con una tolerancia.
!Desde esta subrutina se llamará a las subrutinas dsilus y a dbcg (en ambos casos se envían cachos de estos vectores iw,w).
!---------------------------------------------------------------------------------------------------------------------------------------------
subroutine dslubc (b,n,x,lenw,leniw)	 
use allocatacion
integer*4 crb,cib
parameter (crb=1,cib=11) 
integer*4, dimension(:),allocatable::iw
integer*4 u,ie,ite,itmax,itl,leniw,lenw,n,nt,icol,j,jbgn,jend
integer*4 cdin,cdz,cil,ciu,ciw,cjl,cju,cl,cnc,cnr,cp,cpp,cr,crr,cu,cw,cz,czz,nl,nu
real*8, dimension(:),allocatable::w
real*8 er,tl,b(n),x(n) 
character tiempo*8

allocate(iw(leniw),w(lenw))
	  
!Se muestra por pantalla la hora (con formato hh:mm:ss) cuando se entra aquí
!---------------------------------------------------------------------------
call time(tiempo)
write(6,*) 'Hora:',tiempo

!Inicialización previa de variables:
!-----------------------------------
do u=1,lenw									
w(u)=0.0_8
enddo
do u=1,leniw
iw(u)=0
enddo
!Se toma 'itl=1'. Las posibilidades 'itl=1' y 'itl=2' conllevan a los mismos criterios similares a 'itol=1' y 'itol=2' en la subrutina 
!gradientesbiconjugados. 
!Si 'itl=1', las iteraciones paran cuando la norma 2 (el módulo) del residual dividida por la norma 2 del lado derecho es menor que 
!tol, [Ax-b]/[b]<tol. 
!Si 'itl=2', las iteraciones paran cuando la norma 2 de la inversa de M por el residual divido por la norma 2 de la inversa de M por el 
!término del lado derecho del sistema es menor que tol, [inv(M)*(Ax-b)]/[inv(M)*b]<tol.   
!M es el precondicionador tal que M=matrizsimilar(A). Sin embargo aquí se utilizará la matriz con la diagonal de A. 	  
itl=1
tl=1e-10
itmax=15*n
!nt=ndim-1 es el número de elementos no nulos en A (diagonal y elem no nulos fuera de ella).
nt=cja(n+1)-1
!Se podría usar la magnitud positiva más pequeña para inicializar er: er=tiny(1.0_8).
ite=0
er=0.0_8
ie=0

!Con precondicionamiento (precondicionador LU). 
!Es válido en todos los casos (haya o no ceros en la diagonal).
!Buenos resultados en todos los casos pero parece que en mallas no regulares (backward) SW, NS con pesos BG en general hay inestabilidades 
!del método para altos números de Reynolds donde existiría convergencia del sist no lineal. En este caso (Re alto) funcionaría peor con 
!estabilización y aún peor con Newton. Parece que en mallas regulares con SW, NS no había problema con pesos BG pudiendo incluso resolverse 
!para Re alto cuando el precondicionador diagonal falla (se observó para gran refinamiento) existiendo convergencia del sist. no lineal.   
!Con estabilización podría haber problemas a mayor refinamiento, y con newton con poco refinamiento (Re alto).
!Se tiene un menor número de iteraciones que con método PCGB con precondicionador diagonal.
	 	  
!Comprueba si hay sistema (siempre habrá)
!----------------------------------------
if ((n.lt.1).or.(nt.lt.1)) then
ie=3
write(6,*)'Se tiene el error (buscar por numero):',ie
return
endif

 !Cuenta el número de coeficientes no nulos de la matriz de precondicionamiento ilu.
 !---------------------------------------------------------------------------------- 
 !Cálculo de 'nl,nu'.
  nl=0
  nu=0
  do icol=1,n						 
    !No cuenta en la diagonal.
    jbgn=cja(icol)+1
    jend=cja(icol+1)-1
    if(jbgn.le.jend) then
    !Si hay coeficientes no nulos (diferentes del de la diagonal) en la columna (A está almacenada por columnas).
      do j=jbgn,jend
        if(cia(j).gt.icol) then
		 !Coeficientes de la matriz triangular inferior.
         nl=nl+1
        else
		 !Coeficientes de la matriz triangular superior.
         nu=nu+1
        endif
      enddo
    endif
  enddo

 !A continuación establecen las matrices de trabajo.
 !--------------------------------------------------
 cil=cib
 cjl=cil+n+1
 ciu=cjl+nl
 cju=ciu+nu
 cnr=cju+n+1
 cnc=cnr+n
 ciw=cnc+n

 cl=crb
 cdin=cl+nl
 cu=cdin+n
 cr=cu+nu
 cz=cr+n
 cp=cz+n
 crr=cp+n
 czz=crr+n
 cpp=czz+n
 cdz=cpp+n
 cw=cdz+n

 !Chequea el espacio reservado para los vectores 'w,iw'.
 call dchkw(ciw,leniw,cw,lenw,ie,ite) 
 if (ie.eq.1) then
 write(6,*)'Se tiene el error (buscar por numero):',ie
 return
 endif

 !Si aún no hay error previos se referencia en el vector puntero 'iw' en qué partes de los vectores 'iw,w' irán
 !las matrices que se almacenarán.
 iw(1)=cil
 iw(2)=cjl
 iw(3)=ciu
 iw(4)=cju
 iw(5)=cl
 iw(6)=cdin
 iw(7)=cu
 iw(9)=ciw
 iw(10)=cw

 !A continuación se computa la descomposición lu incompleta.
 !----------------------------------------------------------
 !Sea A, de 10 componentes. A(2:4) es un vector con la primera componente en 2 y la última en 4 de modo que tiene 3 componentes.
 !Cuando se pasa a una subrutina A(2) se envía un trozo del vector A. El vector que se genera en la siguiente tendrá una dimensión 8.
 !Así pues, se envían a la subrutina dsilus trozos de los vectores.
 call dsilus(n,nl,iw(cil),iw(cjl),w(cl),w(cdin),nu,iw(ciu),iw(cju),w(cu),iw(cnr),iw(cnc))	           																                                               
 
 !A continuación se efectúa el algoritmo de los gradientes biconjugados con el precondicionamiento indicado.
 !----------------------------------------------------------------------------------------------------------
 !Se envía a la subrutina dbcg otros trozos de los vectores. Desde ella se llamará de nuevo a varias subrutinas.
 call dbcg(n,b,x,itl,tl,itmax,ite,er,ie,w(cr),w(cz),w(cp),w(crr),w(czz),w(cpp),w(cdz),w,iw,lenw,leniw)
  
 !Se muestra por pantalla el error obtenido en caso de que haya algún problema
 !---------------------------------------------------------------------------- 
 if (ie.eq.0) then
 write(6,*) 'Solucion obtenida en iteracion:',ite
 else
 write(6,*)'Se tiene el error (buscar por numero):',ie
  if (ite.gt.0) then
  write(6,*)'error=',er
  endif
 endif
deallocate(iw,w)
end

!----------------------------------------------------------------------------------------------------------------------------------
!Subrutina DCHKW
!Esta subrutina chequea si la longitud reservada para las matrices de trabajo (iw,w) es suficiente y devuelve el error 
!correspondiente en otro caso. 
!----------------------------------------------------------------------------------------------------------------------------------
subroutine dchkw (lociw,leniw,locw,lenw,ierr,iter) 
integer*4 ierr,iter,leniw,lenw,lociw,locw

!Inicialización previa de variables
!----------------------------------
ierr=0
iter=0
  
!Chequeo de la longitud reservada
!--------------------------------
!Chequea el espacio de la matriz de coeficientes enteros de trabajo.
if(lociw.gt.leniw) then 	 
   ierr=1
   write(6,*) 'Espacio necesario matriz entera:',lociw
   write(6,*) 'Espacio disponible:',leniw
endif
!Chequea el espacio de la matriz de coeficientes reales de trabajo.
if(locw.gt.lenw) then 
   ierr=1
   write(6,*) 'Espacio necesario matriz real:',locw
   write(6,*) 'Espacio disponible:',lenw
endif
end

!--------------------------------------------------------------------------------------------------------------------
!Subrutina DSLUI
!Esta subrutina actúa como una interfaz entre una subrutina y la subrutina que realmente computa x=inv(LDU)b.
!iw(1)=localización del primer coeficiente de 'il' en iw, iw(2)=localización del primer coeficiente de 'jl' en iw.
!iw(3)=localización del primer coeficiente de 'iu' en iw, iw(4)=localización del primer coeficiente de 'ju' en iw.
!iw(5)=localización del primer coeficiente de 'l' en rw,  iw(6)=localización del primer coeficiente de 'dinv' en rw.
!iw(7)=localización del primer coeficiente de 'u' en rw.
!--------------------------------------------------------------------------------------------------------------------
subroutine dslui (n,b,x,rw,iw,lw,liw)
use allocatacion
integer*4 n,lw,liw,nu,nl,iw(liw),locdin,locil,lociu,locjl,locju,locl,locu
real*8 b(n),rw(lw),x(n) 

!Saca la ubicación donde se escribirán las matrices que sostienen la factorización ILU.
!--------------------------------------------------------------------------------------
locil=iw(1)
locjl=iw(2)
lociu=iw(3)
locju=iw(4)
locl=iw(5)
locdin=iw(6)
locu=iw(7)
nl=iw(3)-iw(2)
nu=iw(4)-iw(3)

!Resuelve el sistema LUx=b
!-------------------------
call dslui2(n,b,x,nu,nl,iw(locil),iw(locjl),rw(locl),rw(locdin),iw(lociu),iw(locju),rw(locu))
end

!---------------------------------------------------------------------------------------------------------------------------
!Subrutina DSLUI2
!Subrutina para resolver un sistema de la forma L*D*Ux=b, donde L es una matriz triangular inferior, D es una matriz 
!diagonal y U es una matriz triangular superior.
!'il,jl,l' contendrán la unidad triangular inferior L de la descomposición incompleta de la matriz A, almacenada 
!por filas (formato CSR similar al CSC). 'iu,ju,u' contendrán la unidad triangular superior U de la descomposición 
!incompleta de la matriz A, almacenada por columnas (formato CSC). Esta factorización ILU es computada por la subrutina 
!dsilus. La diagonal (todos sus coeficientes tendrán valor 1) es almacenada.
!--------------------------------------------------------------------------------------------------------------------------- 
subroutine dslui2 (n,b,x,nu,nl,il,jl,l,dinv,iu,ju,u)
integer*4 n,nu,nl,il(n+1),iu(nu),jl(nl),ju(n+1),i,icol,irow,j,jbgn,jend
real*8 b(n),dinv(n),x(n),l(nl),u(nu) 

!Se almacena el término independiente (b) en 'x'
!-----------------------------------------------
do i=1,n
  x(i)=b(i)
enddo

!Resuelve L*y=b, almacenando el resultado (y) en 'x' (L almacenada por filas).
!-----------------------------------------------------------------------------
do irow=2,n
  jbgn=il(irow)
  jend=il(irow+1)-1
  if(jbgn.le.jend) then
    do j=jbgn,jend
      x(irow)=x(irow)-l(j)*x(jl(j))
    enddo
  endif
enddo

!Resuelve D*z=y, almacenando el resultado (z) en 'x'.
!----------------------------------------------------
do i=1,n
  x(i)=x(i)*dinv(i)
enddo

!Resuelve U*x=z, (U almacenada por columnas).
!--------------------------------------------	  
do icol=n,2,-1
  jbgn=ju(icol)
  jend=ju(icol+1)-1
  if(jbgn.le.jend) then
    do j=jbgn,jend
      x(iu(j))=x(iu(j))-u(j)*x(icol)
    enddo
  endif
enddo
end

!------------------------------------------------------------------------------------------------------------------------------
!Subrutina DSILUS
!Subrutina para generar la descomposición incompleta LDU de la matriz siendo L la unidad triangular inferior, U la unidad
!triangular superior. Se almacena la inversa de la diangonal D.	 
!'nl' es el número de coeficientes no nulos en la matriz L, 'nu' es el número de coeficientes no nulos en la matriz U.
!'nrow(i)' es el número de coeficientes no nulos en la fila 'i' de L, 'ncol(i)' es el número de coeficientes no nulos 
!en la columna 'i' de U.
!------------------------------------------------------------------------------------------------------------------------------ 
subroutine dsilus (n,nl,il,jl,l,dinv,nu,iu,ju,u,nrow,ncol)
use allocatacion
integer*4 n,nl,nu,il(n+1),iu(nu),jl(nl),ju(n+1),ncol(n),nrow(n),i,ibgn,icol,iend,indx,indx1,indx2
integer*4 indxc1,indxc2,indxr1,indxr2,irow,itemp,j,jbgn,jend,jtemp,k,kc,kr,kk 
real*8 dinv(n),l(nl),u(nu),temp
      
!El significado de las variables (cuál lleva la posición y cuál lleva la información de donde empieza cada columna) enteras 
!que llevan la matriz "jmatriz y imatriz" cambia entre la matriz triangular inferior y superior.
!'cja,cia,ca' (formato columna) equivale a 'ju,iu,u' (formato columna) y equivale a 'il,jl,l' (formato fila). 

  !Cálculo del número de coeficientes en cada fila de L y en cada columna de U.
  !----------------------------------------------------------------------------
   do i=1,n
     nrow(i)=0
     ncol(i)=0
   enddo
   do icol=1,n
     jbgn=cja(icol)+1
     jend=cja(icol+1)-1
     if(jbgn.le.jend) then
	 !Si hay coeficientes en esa columna.
       do j=jbgn,jend
         if(cia(j).lt.icol) then
          !Cálculo del número de elementos no nulos de la matriz triangular superior en la columna 'icol'.
		  !Se añaden todos para cada 'icol'.
          ncol(icol)=ncol(icol)+1
         else
	      !cia(j) contiene la fila en que está posicionado el coeficiente dentro de la columna 'icol'.
		  !Cálculo del número de elementos no nulos en la columna de la matriz triangular inferior.
	      !Se van añadiendo elementos en la misma fila para varios 'icol'.
          nrow(cia(j))=nrow(cia(j))+1
         endif
       enddo
     endif
   enddo
   
   !Copia la matriz A dentro de las estructuras L y U.
   !--------------------------------------------------
   !Se arreglan los vectores punteros 'ju,il'
   ju(1)=1
   il(1)=1
   do icol=1,n
     il(icol+1)=il(icol)+nrow(icol)
     ju(icol+1)=ju(icol)+ncol(icol)
     nrow(icol)=il(icol)
     ncol(icol)=ju(icol)
   enddo
   !Se escriben los coeficientes y se arreglan los vectores de componentes enteras 'iu,jl'.
   do icol=1,n
     dinv(icol)=ca(cja(icol))
     jbgn=cja(icol)+1
     jend=cja(icol+1)-1
     if(jbgn.le.jend) then
       do j=jbgn,jend
         irow=cia(j)
         if(irow.lt.icol) then
           !Coeficientes del triángulo superior.
           iu(ncol(icol))=irow
           u(ncol(icol))=ca(j)
           ncol(icol)=ncol(icol)+1
         else
           !Coeficientes del triángulo inferior (almacenado por filas).
           jl(nrow(irow))=icol
           l(nrow(irow))=ca(j)
           nrow(irow)=nrow(irow)+1
         endif
       enddo
     endif
   enddo
   !L está almacenada en las primeras 'nl' componentes de 'l,jl,il' y U está almacenada en las primeras 'nu' 
   !componentes de 'u,iu,ju'.

   !Se ordenan las filas de L y las columnas de U (creo que ya están ordenadas y este procedimiento no es necesario)
   !----------------------------------------------------------------------------------------------------------------
   do k=2,n
     jbgn=ju(k)
     jend=ju(k+1)-1
     if(jbgn.lt.jend) then
       do j=jbgn,jend-1
         do i=j+1,jend
           if(iu(j).gt.iu(i)) then
             itemp=iu(j)
             iu(j)=iu(i)
             iu(i)=itemp
             temp=u(j)
             u(j)=u(i)
             u(i)=temp
           endif
         enddo
       enddo
     endif
     ibgn=il(k)
     iend=il(k+1)-1
     if(ibgn.lt.iend) then
       do i=ibgn,iend-1
         do j=i+1,iend
           if(jl(i).gt.jl(j)) then
             jtemp=ju(i)
             ju(i)=ju(j)
             ju(j)=jtemp
             temp=l(i)
             l(i)=l(j)
             l(j)=temp
           endif
         enddo
       enddo
     endif
   enddo

   !Se realiza la descomposición LDU incompleta.
   !--------------------------------------------
   do i=2,n
   !Se evalúa la fila 'i' de L.
     indx1=il(i)
     indx2=il(i+1)-1
     if(indx1.le.indx2) then
     do indx=indx1,indx2
       if(indx.ne.indx1) then
       indxr1=indx1
       indxr2=indx-1
       indxc1=ju(jl(indx))
       indxc2=ju(jl(indx)+1)-1
       if(indxc1.le.indxc2) then
 		  kk=1	
		  kr=jl(indxr1)
		  do while (kk.eq.1)	
		  kk=0
		  kc=iu(indxc1) 
		   if(kr.gt.kc) then
           indxc1=indxc1+1
            if(indxc1.le.indxc2) then     
			kk=1
			endif	      
		   elseif(kr.lt.kc) then
           indxr1=indxr1+1
            if(indxr1.le.indxr2) then
			kr=jl(indxr1)
			kk=1
			endif
           else
           l(indx)=l(indx)-l(indxr1)*dinv(kc)*u(indxc1)
           indxr1=indxr1+1 
           indxc1=indxc1+1 
            if((indxr1.le.indxr2).and.(indxc1.le.indxc2)) then
			kr=jl(indxr1)
			kk=1
			endif
		   endif
 		  enddo
	   endif
	   endif
       l(indx)=l(indx)/dinv(jl(indx))
     enddo
	 endif

     !Se evalúa la columna 'i' de U.
     indx1=ju(i)
     indx2=ju(i+1)-1
     if(indx1.le.indx2) then
     do indx=indx1,indx2
       if(indx.ne.indx1) then
       indxc1=indx1
       indxc2=indx-1	
       indxr1=il(iu(indx))
       indxr2=il(iu(indx)+1)-1
       if(indxr1.le.indxr2) then
		  kk=1
		  kr=jl(indxr1)
		  do while (kk.eq.1)	
		  kk=0
		  kc=iu(indxc1)
		   if(kr.gt.kc) then
           indxc1=indxc1+1
            if(indxc1.le.indxc2) then
			kk=1
			endif
           elseif(kr.lt.kc) then
           indxr1=indxr1+1
            if(indxr1.le.indxr2) then			   
			kr=jl(indxr1)
			kk=1
			endif
           else
           u(indx)=u(indx)-l(indxr1)*dinv(kc)*u(indxc1)
           indxr1=indxr1+1
           indxc1=indxc1+1
            if((indxr1.le.indxr2).and.(indxc1.le.indxc2)) then			   
			kr=jl(indxr1)
			kk=1
			endif
           endif
		  enddo
	   endif
	   endif
       u(indx)=u(indx)/dinv(iu(indx))
     enddo
	 endif

     !Se evalúa el coeficiente 'i' de la diagonal.
     indxr1=il(i)
     indxr2=il(i+1)-1
     if(indxr1.le.indxr2) then  
     indxc1=ju(i)
     indxc2=ju(i+1)-1
     if(indxc1.le.indxc2) then 
	   kk=1
	   kr=jl(indxr1)
	   do while (kk.eq.1)	
	   kk=0
	   kc=iu(indxc1)
		if(kr.gt.kc) then
        indxc1=indxc1+1
         if(indxc1.le.indxc2) then
		 kk=1
		 endif
        elseif(kr.lt.kc) then
        indxr1=indxr1+1
         if(indxr1.le.indxr2) then			   
		 kr=jl(indxr1)
		 kk=1
		 endif
        else
        dinv(i)=dinv(i)-l(indxr1)*dinv(kc)*u(indxc1)
        indxr1=indxr1+1
        indxc1=indxc1+1
         if((indxr1.le.indxr2).and.(indxc1.le.indxc2)) then			
		 kr=jl(indxr1)
		 kk=1
		 endif
        endif
	   enddo  
	 endif
	 endif
   enddo

   !Se reemplazan los elementos diagonales por sus inversos. 
   do i=1,n
    dinv(i)=1.0_8/dinv(i)
   enddo
end

!-------------------------------------------------------------------------------------------------------------------------
!Subrutina DSLUTI
!Esta subrutina actúa como una interfaz entre una subrutina y la subrutina que realmente computa x=inv(transpuesta(LDU))b.
!iw(1)=localización del primer coeficiente de 'il' en 'iw', iw(2)=localización del primer coeficiente de 'jl' en 'iw'.
!iw(3)=localización del primer coeficiente de 'iu' en 'iw', iw(4)=localización del primer coeficiente de 'ju' en 'iw'.
!iw(5)=localización del primer coeficiente de 'l' en 'rw',  iw(6)=localización del primer coeficiente de 'dinv' en 'rw'.
!iw(7)=localización del primer coeficiente de 'u' en 'rw'.
!-------------------------------------------------------------------------------------------------------------------------
subroutine dsluti (n,b,x,rw,iw,lw,liw)
integer*4 n,lw,liw,nu,nl,iw(liw),locdin,locil,lociu,locjl,locju,locl,locu
real*8 b(n),rw(lw),x(n) 

!Saca los punteros a las matrices L,D y U.
!-----------------------------------------
locil=iw(1)
locjl=iw(2)
lociu=iw(3)
locju=iw(4)
locl=iw(5)
locdin=iw(6)
locu=iw(7)
nl=iw(3)-iw(2)
nu=iw(4)-iw(3)

!Resuelve el sistema transpuesta(LDU)x=b
!---------------------------------------
!Los arrays con () dentro de un call son nuevos arrays cuya primera componente es la que va entre paréntesis.  
call dslui4(n,b,x,nu,nl,iw(locil),iw(locjl),rw(locl),rw(locdin),iw(lociu),iw(locju),rw(locu))
end

!-------------------------------------------------------------------------------------------------------------------------
!Subrutina DSLUI4
!Subrutina para resolver un sistema de la forma transpuesta(L*D*U)x=b, donde L es una matriz triangular inferior 
!D es una matriz diagonal y U es una matriz triangular superior.
!-------------------------------------------------------------------------------------------------------------------
subroutine dslui4 (n,b,x,nu,nl,il,jl,l,dinv,iu,ju,u)							   
integer*4 n,nu,nl,il(n+1),iu(nu),jl(nl),ju(n+1),i,icol,irow,j,jbgn,jend
real*8 b(n),dinv(n),x(n),l(nl),u(nu) 

!Se almacena el término independiente (b) en 'x'
!-----------------------------------------------
do i=1,n
  x(i)=b(i)
enddo

!Resuelve transpuesta(U)*y=x, almacenando el resultado (y) en 'x' (U almacenada por columnas).	
!---------------------------------------------------------------------------------------------
do irow=2,n
  jbgn=ju(irow)
  jend=ju(irow+1)-1
  if(jbgn.le.jend) then
    do j=jbgn,jend
      x(irow)=x(irow)-u(j)*x(iu(j))
    enddo
  endif
enddo

!Resuelve D*z=y, almacenando el resultado (z) en 'x'.	
!---------------------------------------------------- 
!'dinv' es la inversa de la matriz diagonal D.
do i=1,n
  x(i)=x(i)*dinv(i)
enddo

!Resuelve transpuesta(L)*x=z (L almacenada por filas).
!-----------------------------------------------------
do icol=n,2,-1
  jbgn=il(icol)
  jend=il(icol+1)-1
  if(jbgn.le.jend) then
    do j=jbgn,jend
      x(jl(j))=x(jl(j))-l(j)*x(icol)
    enddo
  endif
enddo
end

!--------------------------------------------------------------------------------------------------------------------------------
!Subrutina DSMTV
!Subrutina para calcular el producto matriz-vector: y = transpuesta(A)*x. La matriz dispersa A estará almacenada en formato CSC.
!Es equivalente a la subrutina dsprstx (considera formato MSR para A). 
!--------------------------------------------------------------------------------------------------------------------------------
subroutine dsmtv (n,x,y) 	 
use allocatacion
integer*4 n,i,irow
real*8 x(n),y(n) 

!Inicializa a cero el vector resultado.
!--------------------------------------
do i=1,n
  y(i)=0.0_8
enddo

!Multiplica por transpuesta(A).
!------------------------------ 
!Se tiene en cuenta que se tiene transpuesta(A) si se considera que A está almacenada por filas.
do irow=1,n
  do i=cja(irow),cja(irow+1)-1
    y(irow)=y(irow)+ca(i)*x(cia(i))
  enddo
enddo
end

!---------------------------------------------------------------------------------------------------------------------------------
!Subrutina DSMV
!Subrutina para calcular el producto matriz-vector: y = A*x. Es	equivalente a la subrutina dsprsax (considera formato MSR para A).
!---------------------------------------------------------------------------------------------------------------------------------
subroutine dsmv (n,x,y)
use allocatacion
integer*4 n,i,icol
real*8 x(n),y(n) 

!Inicializa a cero el vector resultado.
!--------------------------------------
do i=1,n
  y(i)=0.0_8
enddo

!Multiplica por A.
!-----------------
do icol=1,n
  do i=cja(icol),cja(icol+1)-1
    y(cia(i))=y(cia(i))+ca(i)*x(icol)
  enddo
enddo
end

!---------------------------------------------------------------------------------------------------------------------------
!Subrutina DBCG
!Subrutina para resolver un sistema lineal no simétrico Ax=b usando el método de los gradientes biconjugados precondicionado
!Se llevan a cabo las iteraciones del método. Es equivalente a la subrutina gradientesbiconjugados.
!'w' contiene por este orden: (l,dinv,u,r,z,p,rr,zz,pp,dz) con dimensiones (nl,n,nu,n,n,n,n,n,n,n) respectivamente.
!'iw' contiene por este orden: (índices,il,jl,iu,ju,nrow,ncol) con dimensiones (10,n+1,nl,nu,n+1,n) respectivamente.
!---------------------------------------------------------------------------------------------------------------------------
subroutine dbcg(n,b,x,itol,tol,itmax,iter,err,ierr,r,z,p,rr,zz,pp,dz,rw,iw,lenw,leniw)   
use allocatacion
integer*4 i,k,ierr,iter,itmax,itol,n,lenw,leniw,iw(leniw),isdbcg
real*8 err,tol,b(n),dz(n),p(n),pp(n),r(n),rr(n),rw(lenw),x(n),z(n),zz(n)  
real*8 ak,akden,bk,bkden,bknum,bnrm,tolmin,fuzz   			   

!Inicialización previa de variables
!----------------------------------
iter=0
!Se usa el espacio relativo más pequeño  
fuzz=epsilon(1.0_8)/2.0_8	
tolmin=500.0_8*fuzz
fuzz=fuzz*fuzz
!Chequeo de datos de entrada:
if(tol.lt.tolmin) then
  tol=tolmin
  ierr=4
endif

 !Se calcula el residual inicial y el pseudo-residual 
 !---------------------------------------------------  
 call dsmv(n,x,r)
 do i=1,n
    r(i)=b(i)-r(i)
    rr(i)=r(i)
 enddo
 call dslui(n,r,z,rw,iw,lenw,leniw)
 call dsluti(n,rr,zz,rw,iw,lenw,leniw)

 !Si se cumple el criterio de parada se sale de la subrutina sin hacer iteración alguna.
 if(isdbcg(n,b,itol,tol,iter,err,ierr,r,z,dz,rw,iw,bnrm,lenw,leniw).ne.0) then
 return
 endif
 if(ierr.ne.0) then
 return
 endif

 !Comienzo de las iteraciones del método
 !--------------------------------------
 do k=1,itmax
   iter=k

   !Calcula el coeficiente 'bk' y los vectores de dirección 'p' y 'pp'.
   bknum=0.0_8
   !Se calcula el producto interior de dos vectores.
   bknum=dot_product(z,rr)
   if(abs(bknum).le.fuzz) then
     ierr=5
     return
   endif
   if(iter.eq.1) then
     !Se hace la copia un vector.
     do i=1,n
     p(i)=z(i)
     pp(i)=zz(i)
     enddo
   else
     bk=bknum/bkden
     do i=1,n
       p(i)=z(i)+bk*p(i)
       pp(i)=zz(i)+bk*pp(i)
     enddo
   endif
   bkden=bknum

   !Calcula el coeficiente 'ak', nueva iteración de 'x', nuevos residuales 'r' y 'rr' y 
   !nuevos pseudo-residuales 'z' y 'zz'.
   call dsmv(n,p,z)
   akden=0.0_8
   !Se calcula el producto interior de dos vectores.
   akden=dot_product(pp,z)	 
   ak=bknum/akden
   if(abs(akden).le.fuzz) then
     ierr=6
     return
   endif
   !Se calcula la suma de una constante 'ak' multiplicada por un vector y un vector.
   if (ak.ne.0.0_8) then
    do i=1,n
    x(i)=x(i)+ak*p(i)
    r(i)=r(i)-ak*z(i)
    enddo
   endif
   call dsmtv(n,pp,zz)
   !Se calcula la suma de una constante (ak) multiplicada por un vector y un vector
   if (ak.ne.0.0_8) then
    do i=1,n
    rr(i)=rr(i)-ak*zz(i)
    enddo
   endif
   call dslui(n,r,z,rw,iw,lenw,leniw)
   call dsluti(n,rr,zz,rw,iw,lenw,leniw)

   !Evalúa el criterio de parada.
   if (isdbcg(n,b,itol,tol,iter,err,ierr,r,z,dz,rw,iw,bnrm,lenw,leniw).ne.0) then 
   return
   endif
 enddo

!Criterio de parada no satisfecho.
iter=itmax+1
ierr=2
end   

!-------------------------------------------------------------------------------------------------------------------------------------
!Función ISDBCG
!Esta subrutina evalúa la parada para el esquema iterativo BCG. Devuelve un valor no nulo si el error estimado (tipo de 
!estimación determinada por 'itol') es menor que la tolerancia especificada (tol).
!'bnrm' es la norma del término del lado derecho del sistema considerado (depende de 'itol'). Es calculado sólo en la primera llamada. 
!Valores que devuelve la función:
! 0 : El estimador del error no es menor que la tolerancia especificada. Las iteraciones deben continuar.
! 1 : El estimador del error es menor que la tolerancia especificada. El proceso iterativo se considera completado.
!-------------------------------------------------------------------------------------------------------------------------------------
function isdbcg(n,b,itol,tol,iter,err,ierr,r,z,dz,rw,iw,bnrm,lenw,leniw)
use allocatacion
integer*4 ierr,iter,itol,n,lenw,leniw,iw(leniw),isdbcg
real*8 bnrm,err,tol,b(n),dz(n),r(n),rw(lenw),z(n),dnrm2

302 format (A7,I8,5X,A7,E11.4E2)

 !Cálculo del error estimado
 !--------------------------
 isdbcg=0
 if (itol.eq.1) then
   !err = ||residual||/||b|| (el símbolo ||^|| indica la norma 2 de ^).
   if (iter.eq.0) then
   bnrm=dnrm2(n,b)
   endif
   err=dnrm2(n,r)/bnrm  
 elseif (itol.eq.2) then
   !err = ||inv(M)*residual||/||inv(M)*b||.
   if (iter.eq.0) then
     call dslui(n,b,dz,rw,iw,lenw,leniw)
     bnrm=dnrm2(n,dz)  
   endif
   err=dnrm2(n,z)/bnrm	
 else
   !Si se entra aquí, itol no tiene un valor correcto.
   ierr=3
 endif

 !Escritura por pantalla de iteracion y del error.
 !------------------------------------------------
 if ((iter.eq.100).or.(iter.eq.1000).or.(iter.eq.n).or.(iter.eq.7*n).or.(iter.eq.15*n)) then
 write(6,302) ' iter =',iter,'error =',err
 endif

 !Se evalúa si el error es mayor o menor que la tolerancia.
 !---------------------------------------------------------
 if (err.le.tol) then
 isdbcg=1
 endif
end

!--------------------------------------------------------------------------------------------------------------------------------------------
!Función DNRM2
!Esta función calcula la longitud Euclidiana (norma L2) de un vector (la raiz del cuadrado de sus componentes).
!'n' es el número de elementos en el vector, 'dx' es el array donde está el vector y tiene dimensión 'n' (los coeficientes no están 
!colocados espaciados). Es equivalente a la función snrm.
!Se trabaja con dos contantes 'cutlo' y 'cuthi', y se evalúan cuatro fases.
!La fase 1 escanea componentes nulas. Se va a la fase 2 cuando una componente (distinta de cero) es menor o igual que 'cutlo'.
!Se va a la fase 3 cuando una componente es mayor que 'cutlo'. Se va a la fase 4 cuando una componente es mayor o igual que 'cuthi/n'.
!--------------------------------------------------------------------------------------------------------------------------------------------
function dnrm2 (n,dx)
integer*4 n,i,j  
real*8 dx(n),cutlo,cuthi,hitest,sum,xmax,dnrm2 
  
!Inicialización previa de variables
!----------------------------------   
cutlo=sqrt(tiny(1.0_8)/epsilon(1.0_8))
cuthi=sqrt(huge(1.0_8))

if (n.le.0) then
  dnrm2=0.0_8 
  return
endif         			
sum=0.0_8 
	   
   !Comienzo del bucle principal
   !----------------------------
   i=1
   xmax=0.0_8
   !Fase 1. 'sum' es cero
   do j=i,n			   
    if (dx(j).eq.0.0_8) then	    
     if (j.eq.n) then   			   
	 dnrm2=xmax*sqrt(sum)
     !Sale de la subrutina con dnrm2=0 (no hay valores no nulos en el vector).
     return	
     endif
    else
    !Elimina los primeros ceros.	                            
    i=j
    exit
    endif
   enddo
	  
   !Evalúa para la primera componente no nula.
   if (abs(dx(i)).gt.cutlo) then	      
   !Establece hitest=cuthi/n.	  
   hitest=cuthi/n
   !Fase 3. 'sum' es de rango medio. Sin escala.
   do j=i,n	                         
    if (abs(dx(j)).ge.hitest) then
    !Se prepara para la fase 4.
    i=j
    sum=(sum/dx(i))/dx(i)
    xmax=abs(dx(i))
    sum=sum+(dx(i)/xmax)**2.0_8
    goto 200
    endif
    sum=sum+dx(j)**2.0_8
   enddo
   dnrm2=sqrt(sum)
   !Sale de la subrutina.
   return                             
   endif

   !Se prepara para la fase 2.
   !Primera componente no nula y menor que cutlo (valor muy pequeño).
   xmax=abs(dx(i))			         
   sum=sum+(dx(i)/xmax)**2.0_8
   i=i+1		                         
   if (i.le.n) then	                 
   !Fase 2. 'sum' es pequeño. Se escala para evitar underflow destructivo. 
70  if (abs(dx(i)).gt.cutlo) then     
    !Se prepara para la fase 3.
	sum=(sum*xmax)*xmax
    !Se establece hitest=cuthi/n                                         
    hitest=cuthi/n
    !Fase 3. 'sum' es de rango medio. Sin escala.
    do j=i,n                          
     if (abs(dx(j)).ge.hitest) then	   
     !Se prepara para la fase 4.
	 i=j
     sum=(sum/dx(i))/dx(i)
     xmax=abs(dx(i))
     sum=sum+(dx(i)/xmax)**2.0_8
     goto 200
     endif
    sum=sum+dx(j)**2.0_8
    enddo
    dnrm2=sqrt(sum)
    return
    endif
    !Código común para fases 2 y 4. En fase 4 'sum' es grande. Se escala para evitar overflow.
    if (abs(dx(i)).le.xmax) then
    sum=sum+(dx(i)/xmax)**2.0_8
    i=i+1	                         
     if (i.le.n) then                 
     goto 70
	 else
     !Fin del buble principal. Cálculo de la raíz cuadrada y ajuste de la escala.		
 	 dnrm2=xmax*sqrt(sum)
	 return
	 endif
    endif
    sum=1.0_8+sum*(xmax/dx(i))**2.0_8 
    xmax=abs(dx(i))
    i=i+1                             
     if (i.le.n) then	             
	 goto 70 
     else
     !Fin del buble principal. Cálculo de la raíz cuadrada y ajuste de la escala.
     dnrm2=xmax*sqrt(sum)
     return
     endif
   
   else
   !Fin del buble principal. Cálculo de la raíz cuadrada y ajuste de la escala.
   dnrm2=xmax*sqrt(sum)
   return
   endif

200 do j=i+1,n                         
    !Código común para fases 2 y 4. En fase 4 sum es grande. Se escala para evitar overflow. 	   
	 if (abs(dx(j)).le.xmax) then	 
      sum=sum+(dx(j)/xmax)**2.0_8	 
      cycle                            
     endif
     sum=1.0_8+sum*(xmax/dx(j))**2.0_8	
     xmax=abs(dx(j))					                                
    enddo                                   
	     
!Fin del buble principal. Cálculo de la raíz cuadrada y ajuste de la escala.
dnrm2=xmax*sqrt(sum)
end

!----------------------------------------------------------------------------------------------------- 
!FIN DE PROGRAMA
!-----------------------------------------------------------------------------------------------------