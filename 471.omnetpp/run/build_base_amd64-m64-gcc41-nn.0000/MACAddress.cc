/*
 * Copyright (C) 2003 CTIE, Monash University
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

#ifdef _MSC_VER
#pragma warning(disable:4786)
#endif

#if defined(SPEC_CPU)
#include <string.h>
#endif

#include <ctype.h>
#include "MACAddress.h"




//
// Converts hex string into a byte array 'destbuf'. Destbuf is 'size'
// chars long -- if hext string is shorter, destbuf is filled with zeros;
// if destbuf is longer, it is truncated. Non-hex characters are
// discarded before conversion. Returns number of bytes converted from hex.
//
static int hextobin(const char *hexstr, unsigned char *destbuf, int size)
{
    int k=0;
    const char *s = hexstr;
    for (int pos=0; pos<size; pos++)
    {
        if (!s || !*s)
        {
            destbuf[pos] = 0;
        }
        else
        {
            while (*s && !isxdigit(*s)) s++;
            if (!*s) {destbuf[pos]=0; continue;}
            unsigned char d = isdigit(*s) ? (*s-'0') : islower(*s) ? (*s-'a'+10) : (*s-'A'+10);
            d = d<<4;
            s++;

            while (*s && !isxdigit(*s)) s++;
            if (!*s) {destbuf[pos]=0; continue;}
            d += isdigit(*s) ? (*s-'0') : islower(*s) ? (*s-'a'+10) : (*s-'A'+10);
            s++;

            destbuf[pos] = d;
            k++;
        }
    }
    return k;
}


MACAddress::MACAddress() : MACAddress_Base()
{
    address[0]=address[1]=address[2]=address[3]=address[4]=address[5]=0;
}

MACAddress::MACAddress(const char *hexstr) : MACAddress_Base()
{
    setAddress(hexstr);
}

MACAddress& MACAddress::operator=(const MACAddress& other)
{
    MACAddress_Base::operator=(other);
    memcpy(address, other.address, MAC_ADDRESS_BYTES);
    return *this;
}

unsigned int MACAddress::getAddressArraySize() const
{
    return 6;
}

unsigned char MACAddress::getAddress(unsigned int k) const
{
    if (k>=6) throw new cException("Array of size 6 indexed with %d", k);
    return address[k];
}

void MACAddress::setAddress(unsigned int k, unsigned char addrbyte)
{
    if (k>=6) throw new cException("Array of size 6 indexed with %d", k);
    address[k] = addrbyte;
}

void MACAddress::setAddress(const char *hexstr)
{
    if (!hexstr)
        throw new cException("MACAddress::setAddress(const char *): got null pointer");
    if (hextobin(hexstr, address, MAC_ADDRESS_BYTES)!=MAC_ADDRESS_BYTES)
        throw new cException("MACAddress::setAddress(const char *): hex string \"%s\" too short, should be 12 hex digits", hexstr);
}

void MACAddress::setAddressBytes(unsigned char *addrbytes)
{
    memcpy(address, addrbytes, MAC_ADDRESS_BYTES);
}

void MACAddress::setBroadcast()
{
    address[0]=address[1]=address[2]=address[3]=address[4]=address[5]=0xff;
}

bool MACAddress::isBroadcast() const
{
    return address[0]==0xff && address[1]==0xff && address[2]==0xff && address[3]==0xff &&
           address[4]==0xff && address[5]==0xff;
}

bool MACAddress::isEmpty() const
{
    return !(address[0] || address[1] || address[2] || address[3] || address[4] || address[5]);
}

const char *MACAddress::toHexString(char *buf) const
{
    char *s = buf;
    for (int i=0; i<MAC_ADDRESS_BYTES; i++, s+=3)
        sprintf(s,"%2.2X:",address[i]);
    *(s-1)='\0';
    return buf;
}

bool MACAddress::equals(const MACAddress& other) const
{
    return memcmp(address, other.address, MAC_ADDRESS_BYTES)==0;
}

int MACAddress::compareTo(const MACAddress& other) const
{
    return memcmp(address, other.address, MAC_ADDRESS_BYTES);
}


