use super::*;

#[test]
fn the_digests_of_the_rfc_test_suite() {
    for (text, want) in [
        ("", "d41d8cd98f00b204e9800998ecf8427e"),
        ("a", "0cc175b9c0f1b6a831c399e269772661"),
        ("abc", "900150983cd24fb0d6963f7d28e17f72"),
        ("message digest", "f96b697d7cb7938d525a2f31aaf161d0"),
        (
            "abcdefghijklmnopqrstuvwxyz",
            "c3fcd3d76192e4007dfb496cca67e13b",
        ),
        (
            "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789",
            "d174ab98d277d9f5a5611c2c9f419d9f",
        ),
        (
            "12345678901234567890123456789012345678901234567890123456789012345678901234567890",
            "57edf4a22be3c955ac49da2e2107b67a",
        ),
    ] {
        assert_eq!(hex(text.as_bytes()), want, "{text}");
    }
}

#[test]
fn lengths_around_a_block_are_padded_right() {
    // 55, 56, 57, 63, 64 and 65 bytes straddle where the padding of 1 + length needs another block
    let a: Vec<String> = [55, 56, 57, 63, 64, 65, 119, 120, 128]
        .iter()
        .map(|&n| hex(&vec![b'a'; n]))
        .collect();
    assert_eq!(
        a.iter().collect::<std::collections::HashSet<_>>().len(),
        a.len()
    );
    assert_eq!(hex(&[b'a'; 55]), "ef1772b6dff9a122358552954ad0df65");
    assert_eq!(hex(&[b'a'; 56]), "3b0c8ac703f828b04c6c197006d17218");
    assert_eq!(hex(&[b'a'; 64]), "014842d480b571495a4a0363793f7367");
    assert_eq!(hex(&vec![b'a'; 1000]).len(), 32);
}
