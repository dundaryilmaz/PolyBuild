
  �+  e   k820309    a          16.0        T�[                                                                                                           
       mcmodule.f90 MCMODULE                                                     
                                                           
                                                           
                         @               �                '�                   #N_BEADS    #CHAIN_INDEX    #N_MERS    #MER_ID    #CLOSE_FLAG 	   #MER_FLAG 
   #HEADH_ATOM_INDEX    #TAILH_ATOM_INDEX    #FUNCTIONAL    #HTYPE    #CHAINTYPE    #COORD    #HEAD    #TAIL    #MERS    #BONDLIST .   #PUT_MONOMER /   #PUT_ALL_MONOMERS 3   #CHECK_TWO_MERS 6   #PRINT_BEADS_ONLY ;   #REARRANGE_MER A   #PRINT_POLYMER_CHAIN G   #MAP_ATYPES L   #TRANSLATE_CHAIN O                �                                                               �                                                              �                                                              �                                                            �                               	                                         &                                                      �                               
            X                             &                                                        �                                    �                          �                                    �                          �                                    �       	                   �                                    �       
                   �                                           �                         �                                          �                 
            &                   &                                                        �                                                          
  p          p            p                                       �                                          0                
  p          p            p                                     �                                           H                  #MONOMER              &                                                             @              @                '                   #MONOMER_NAME    #DELETED    #ATOM_INDEX    #N_MER_ATOMS    #HEADC    #TAILC    #ATYPE    #SYMBOL    #COORD    #MASS    #CHARGE    #BASE_VECS     #MER_LENGTH !   #MER_WIDTH "   #HEADH #   #TAILH $   #HEADCHARGE %   #TAILCHARGE &   #HTYPE '   #RADICAL (   #BONDLIST )   #CREATE *                �                                                                    �                                                                        &                                                      �                                           X                             &                                                        �                                    �                          �                                    �                          �                                    �                        �                                           �                             &                                           .           �                                          �                             &                                                              �                                          @             	   
            &                   &                                                      �                                          �             
   
            &                                                      �                                          �                
            &                                                        �                                           0                
  p          p          p            p          p                                       �                              !     `         
                �                              "     h         
                �                              #            p                
  p          p            p                                       �                              $            �                
  p          p            p                                       �                              %     �         
                �                              &     �         
                �                               '     �                         �                               (     �                       �                               )            �                            &                   &                                           1         �   � $                     �      *                  #READ_MONOMER +   %         @                                 +                           #THIS ,   #FILENAME -                                             ,                   #MONOMER                                              -                     1            �                               .            �                            &                   &                                           1         �   � $                     �      /                  #PUT_MONOMER 0   %         @                                 0                           #THIS 1   #NTH 2                                             1     �              #CHAIN                                               2            1         �   � $                     �      3                  #PUT_ALL_MONOMERS 4   %         @                                 4                           #THIS 5                                             5     �              #CHAIN    1         �   � $                     �      6                  #CHECK_TWO_MERS 7   %         @                                 7                    
       #THIS 8   #A 9   #B :                                             8     �              #CHAIN                                               9                                                       :            1         �   � $                     �      ;                  #PRINT_BEADS_ONLY <   %         @                                 <                           #THIS =   #FILENAME >   #CHAIN_POS ?   #INDX @                                             =     �              #CHAIN                                              >                     1                                           ?                     1                                            @            1         �   � $                     �      A                  #REARRANGE_MER B   %         @                                 B                    
       #THIS C   #NTH_MER D   #DIST1 E   #DIST2 F                                             C     �              #CHAIN                                               D                                                      E     
                                                 F     
       1         �   � $                     �      G                  #PRINT_POLYMER_CHAIN H   %         @                                 H                           #THIS I   #FILENAME J   #INDX K                                             I     �              #CHAIN                                               J                     1                                            K            1         �   � $                     �      L                  #MAP_ATYPES M   %         @                                 M                           #THIS N                                             N     �              #CHAIN    1         �   � $                     �      O                  #TRANSLATE_CHAIN P   %         @                                 P                           #THIS Q   #VEC R                                             Q     �              #CHAIN                                              R                   
     p          p            p                                                                     S                                   &                                           #MONOMER    %         @                                T                    
       #A U   #B V                                             U                   
 
    p          p            p                                                                    V                   
     p          p            p                          (         `                                W                                    
    p          p            p                          %         @                                X                           #X Y                                             Y                   
     p          p            p                          (         `                                Z                                   
    #X [   p          p            p                                                                    [                   
     p          p            p                          (         `                                \                   	                
    #X ]   p          p            p                                                                    ]                   
     p          p            p                          %         @                                 ^                    
       #CHAINA _             D @                               _     �              #CHAIN    %         @                                 `                    
       #CHAINA a   #CHAINB b             D @                               a     �              #CHAIN              D @                               b     �              #CHAIN    %         @                                 c                           #MC_CHAIN d             D @                               d     �              #CHAIN       �         fn#fn !   �   @   J   GLOBAL_VARIABLES    �   @   J   CHECKBOX    >  @   J   POLYMER_TYPE #   ~  �      CHAIN+POLYMER_TYPE +   Q  H   a   CHAIN%N_BEADS+POLYMER_TYPE /   �  H   a   CHAIN%CHAIN_INDEX+POLYMER_TYPE *   �  H   a   CHAIN%N_MERS+POLYMER_TYPE *   )  H   a   CHAIN%MER_ID+POLYMER_TYPE .   q  �   a   CHAIN%CLOSE_FLAG+POLYMER_TYPE ,     �   a   CHAIN%MER_FLAG+POLYMER_TYPE 4   �  H   a   CHAIN%HEADH_ATOM_INDEX+POLYMER_TYPE 4   �  H   a   CHAIN%TAILH_ATOM_INDEX+POLYMER_TYPE .   )  H   a   CHAIN%FUNCTIONAL+POLYMER_TYPE )   q  H   a   CHAIN%HTYPE+POLYMER_TYPE -   �  P   a   CHAIN%CHAINTYPE+POLYMER_TYPE )   	  �   a   CHAIN%COORD+POLYMER_TYPE (   �  �   a   CHAIN%HEAD+POLYMER_TYPE (   Q  �   a   CHAIN%TAIL+POLYMER_TYPE (   �  �   a   CHAIN%MERS+POLYMER_TYPE %   �	  t      MONOMER+POLYMER_TYPE 2     P   a   MONOMER%MONOMER_NAME+POLYMER_TYPE -   R  �   a   MONOMER%DELETED+POLYMER_TYPE 0   �  �   a   MONOMER%ATOM_INDEX+POLYMER_TYPE 1   z  H   a   MONOMER%N_MER_ATOMS+POLYMER_TYPE +   �  H   a   MONOMER%HEADC+POLYMER_TYPE +   
  H   a   MONOMER%TAILC+POLYMER_TYPE +   R  �   a   MONOMER%ATYPE+POLYMER_TYPE ,   �  �   a   MONOMER%SYMBOL+POLYMER_TYPE +   �  �   a   MONOMER%COORD+POLYMER_TYPE *   .  �   a   MONOMER%MASS+POLYMER_TYPE ,   �  �   a   MONOMER%CHARGE+POLYMER_TYPE /   V  �   a   MONOMER%BASE_VECS+POLYMER_TYPE 0     H   a   MONOMER%MER_LENGTH+POLYMER_TYPE /   Z  H   a   MONOMER%MER_WIDTH+POLYMER_TYPE +   �  �   a   MONOMER%HEADH+POLYMER_TYPE +   >  �   a   MONOMER%TAILH+POLYMER_TYPE 0   �  H   a   MONOMER%HEADCHARGE+POLYMER_TYPE 0   "  H   a   MONOMER%TAILCHARGE+POLYMER_TYPE +   j  H   a   MONOMER%HTYPE+POLYMER_TYPE -   �  H   a   MONOMER%RADICAL+POLYMER_TYPE .   �  �   a   MONOMER%BONDLIST+POLYMER_TYPE ,   �  Z   a   MONOMER%CREATE+POLYMER_TYPE *      h       READ_MONOMER+POLYMER_TYPE /   h  U   a   READ_MONOMER%THIS+POLYMER_TYPE 3   �  L   a   READ_MONOMER%FILENAME+POLYMER_TYPE ,   	  �   a   CHAIN%BONDLIST+POLYMER_TYPE /   �  Y   a   CHAIN%PUT_MONOMER+POLYMER_TYPE )     c       PUT_MONOMER+POLYMER_TYPE .   q  S   a   PUT_MONOMER%THIS+POLYMER_TYPE -   �  @   a   PUT_MONOMER%NTH+POLYMER_TYPE 4     ^   a   CHAIN%PUT_ALL_MONOMERS+POLYMER_TYPE .   b  Z       PUT_ALL_MONOMERS+POLYMER_TYPE 3   �  S   a   PUT_ALL_MONOMERS%THIS+POLYMER_TYPE 2     \   a   CHAIN%CHECK_TWO_MERS+POLYMER_TYPE ,   k  h       CHECK_TWO_MERS+POLYMER_TYPE 1   �  S   a   CHECK_TWO_MERS%THIS+POLYMER_TYPE .   &  @   a   CHECK_TWO_MERS%A+POLYMER_TYPE .   f  @   a   CHECK_TWO_MERS%B+POLYMER_TYPE 4   �  ^   a   CHAIN%PRINT_BEADS_ONLY+POLYMER_TYPE .     �       PRINT_BEADS_ONLY+POLYMER_TYPE 3   �  S   a   PRINT_BEADS_ONLY%THIS+POLYMER_TYPE 7   �  L   a   PRINT_BEADS_ONLY%FILENAME+POLYMER_TYPE 8   $  L   a   PRINT_BEADS_ONLY%CHAIN_POS+POLYMER_TYPE 3   p  @   a   PRINT_BEADS_ONLY%INDX+POLYMER_TYPE 1   �  [   a   CHAIN%REARRANGE_MER+POLYMER_TYPE +     }       REARRANGE_MER+POLYMER_TYPE 0   �  S   a   REARRANGE_MER%THIS+POLYMER_TYPE 3   �  @   a   REARRANGE_MER%NTH_MER+POLYMER_TYPE 1     @   a   REARRANGE_MER%DIST1+POLYMER_TYPE 1   [  @   a   REARRANGE_MER%DIST2+POLYMER_TYPE 7   �  a   a   CHAIN%PRINT_POLYMER_CHAIN+POLYMER_TYPE 1   �  r       PRINT_POLYMER_CHAIN+POLYMER_TYPE 6   n  S   a   PRINT_POLYMER_CHAIN%THIS+POLYMER_TYPE :   �  L   a   PRINT_POLYMER_CHAIN%FILENAME+POLYMER_TYPE 6      @   a   PRINT_POLYMER_CHAIN%INDX+POLYMER_TYPE .   M   X   a   CHAIN%MAP_ATYPES+POLYMER_TYPE (   �   Z       MAP_ATYPES+POLYMER_TYPE -   �   S   a   MAP_ATYPES%THIS+POLYMER_TYPE 3   R!  ]   a   CHAIN%TRANSLATE_CHAIN+POLYMER_TYPE -   �!  c       TRANSLATE_CHAIN+POLYMER_TYPE 2   "  S   a   TRANSLATE_CHAIN%THIS+POLYMER_TYPE 1   e"  �   a   TRANSLATE_CHAIN%VEC+POLYMER_TYPE '   �"  �       BASE_MERS+POLYMER_TYPE "   �#  ^       PBC_DIST+CHECKBOX $   �#  �   a   PBC_DIST%A+CHECKBOX $   �$  �   a   PBC_DIST%B+CHECKBOX $   %  �       RANDOM_UNIT_VEC+VEC *   �%  W       INSIDEVOIDS_BUFF+CHECKBOX ,   &  �   a   INSIDEVOIDS_BUFF%X+CHECKBOX    �&  �       PBC_X+CHECKBOX !   R'  �   a   PBC_X%X+CHECKBOX $   �'  �       RETURN_BOX+CHECKBOX &   �(  �   a   RETURN_BOX%X+CHECKBOX '   %)  \       SELF_CHAIN_INTERACTION .   �)  S   a   SELF_CHAIN_INTERACTION%CHAINA &   �)  h       TWO_CHAIN_INTERACTION -   <*  S   a   TWO_CHAIN_INTERACTION%CHAINA -   �*  S   a   TWO_CHAIN_INTERACTION%CHAINB    �*  ^       CHAIN_MOVE $   @+  S   a   CHAIN_MOVE%MC_CHAIN 